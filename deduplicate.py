#!/usr/bin/env python

import argparse
import gzip
import os
import sys
import zfile

def deduplicate(
	inputPaths,
	gprobsFormat,
	doseFormat,
	bglPaths,
	outputPrefix,
	dupePrefix,
	compressed,
	outCompressed,
	force,
	verbose
):
	files = range(len(inputPaths))
	
	print "opening input files ..."
	inputFiles = list()
	inputDupes = list()
	inputCols = list()
	bglFiles = list()
	outputPaths = list()
	dupePaths = list()
	for f in files:
		# find and open input file
		if (compressed != None) and (os.path.exists(inputPaths[f] + compressed)):
			inputFiles.append(gzip.open(inputPaths[f] + compressed,'rb'))
		elif os.path.exists(inputPaths[f]):
			inputFiles.append(open(inputPaths[f],'rU'))
		else:
			exit("ERROR: could not find input file #%d (%s)" % (f,inputPaths[f]))
		
		# initialize per-file data
		inputDupes.append(set())
		inputCols.append(0)
		
		# find and open .bgl file
		if (compressed != None) and (os.path.exists(bglPaths[f] + compressed)):
			bglFiles.append(gzip.open(bglPaths[f] + compressed,'rb'))
		elif os.path.exists(bglPaths[f]):
			bglFiles.append(open(bglPaths[f],'rU'))			
		else:
			exit("ERROR: could not find .bgl file #%d (%s)" % (f,bglPaths[f]))
		
		# pick and check output path
		if os.path.isdir(outputPrefix):
			outputPaths.append(os.path.join(outputPrefix, os.path.basename(inputPaths[f])))
		else:
			outputPaths.append(outputPrefix + os.path.basename(inputPaths[f]))
		if compressed:
			outputPaths[-1] += compressed
		if (not force) and os.path.exists(outputPaths[f]):
			exit("ERROR: deduplicated output file #%d (%s) already exists" % (f,outputPaths[f]))
		
		# pick and check dupe path
		if os.path.isdir(dupePrefix):
			dupePaths.append(os.path.join(dupePrefix, os.path.basename(inputPaths[f])))
		else:
			dupePaths.append(dupePrefix + os.path.basename(inputPaths[f]))
		if compressed:
			dupePaths[-1] += compressed
		if (not force) and os.path.exists(dupePaths[f]):
			exit("ERROR: duplicate output file #%d (%s) already exists" % (f,dupePaths[f]))
	#foreach file
	print "... OK"
	
	print "comparing headers ..."
	sampleFirst = {}
	sampleDupes = []
	outputLine = []
	dupeLine = []
	for f in files:
		# read and deduplicate header lines
		outputLine.append("marker alleleA alleleB")
		dupeLine.append("marker alleleA alleleB")
		inputCols[f] = len(inputFiles[f].next().strip().split())
		words = bglFiles[f].next().strip().split()
		w = 2
		while w < len(words):
			s = int((w / 2) - 1)
			sample = words[w]
			if sample in sampleFirst:
				sampleDupes.append( sampleFirst[sample]+(f,s) )
				inputDupes[f].add(s)
				h = " %s(%d/%d)" % (sample,sampleFirst[sample][0],f)
				dupeLine[f] += (h*3 if gprobsFormat else h)
			else:
				sampleFirst[sample] = (f,s)
				h = " %s" % (sample,)
				outputLine[f] += (h*3 if gprobsFormat else h)
			w += 2
		outputLine[f] += "\n"
		dupeLine[f] += "\n"
		
		# verify column counts
		c = 3 + ((s + 1) * (3 if gprobsFormat else 1))
		if inputCols[f] != c:
			exit("ERROR: input file #%d (%s) has %d columns, expected %d" % (f,inputPaths[f],inputCols[f],c))
		
		bglFiles[f].close()
		print "  file #%d (%s): %d samples (%d new, %d duplicates)" % (f,inputPaths[f],s,s-len(inputDupes[f]),len(inputDupes[f]))
	#foreach file
	print "... OK"
	
	outputFiles = list()
	dupeFiles = list()
	markerA1 = dict()
	markerA2 = dict()
	for f in files:
		print "deduplicating file #%d (%s) ..." % (f,inputPaths[f])
		
		# open output file
		if outCompressed:
			outputFiles.append(gzip.open(outputPaths[f], 'wb', compresslevel=6))
		else:
			outputFiles.append(open(outputPaths[f], 'w'))
		outputFiles[f].write(outputLine[f])
		
		# open dupe file, if needed
		if inputDupes[f]:
			if outCompressed:
				dupeFiles.append(gzip.open(dupePaths[f], 'wb', compresslevel=6))
			else:
				dupeFiles.append(open(dupePaths[f], 'w'))
			dupeFiles[f].write(dupeLine[f])
		else:
			dupeFiles.append(None)
		
		# process lines
		line = 1
		try:
			while True:
				line += 1
				words = inputFiles[f].next().rstrip("\r\n").split()
				if len(words) != inputCols[f]:
					exit("ERROR: expected %d input columns in file #%d, found %d at line %d" % (inputCols[f],f,len(words),line))
				# for the first file, store the marker's allele order and write straight through
				if f == 0:
					marker = words[0]
					if marker in markerA1:
						exit("ERROR: duplicate marker '%s' on line %d" % (marker,line))
					markerA1[marker] = words[1]
					markerA2[marker] = words[2]
				else:
					# for repeat markers in later files, compare allele order and swap if necessary
					marker = words[0]
					if marker in markerA1:
						if markerA1[marker] == words[2] and markerA2[marker] == words[1]:
							words[1],words[2] = words[2],words[1]
							if gprobsFormat:
								for w in xrange(3,len(words),3):
									words[w],words[w+2] = words[w+2],words[w]
							else:
								for w in xrange(3,len(words)):
									words[w] = "%f" % (2.0 - float(words[w]))
							print "  WARNING: swapped allele order for marker '%s' on line %d" % (marker,line)
						elif markerA1[marker] != words[1] or markerA2[marker] != words[2]:
							exit("ERROR: allele mismatch for marker '%s' on line %d: expected %s/%s, found %s/%s" % (marker,line,markerA1[marker],markerA2[marker],words[1],words[2]))
					else:
						markerA1[marker] = words[1]
						markerA2[marker] = words[2]
				outputFiles[f].write(" ".join(words[w] for w in xrange(0,len(words)) if int((w-3)/(3 if gprobsFormat else 1)) not in inputDupes[f]))
				outputFiles[f].write("\n")
				if inputDupes[f]:
					dupeFiles[f].write(" ".join(words[w] for w in xrange(0,len(words)) if (w < 3 or int((w-3)/(3 if gprobsFormat else 1)) in inputDupes[f])))
					dupeFiles[f].write("\n")
			#forever, until exception
		except StopIteration:
			pass
		
		outputFiles[f].close()
		if inputDupes[f]:
			dupeFiles[f].close()
		print "... OK: %d lines" % line
	#foreach file
#deduplicate()


if __name__ == "__main__":
	versName,versMaj,versMin,versRev,versDate = sys.argv[0],0,1,0,'2012-05-21'
	version = "%d.%d.%d (%s)" % (versMaj, versMin, versRev, versDate)
	
	# define usage
	parser = argparse.ArgumentParser(
		#formatter_class=argparse.RawDescriptionHelpFormatter,
		description="%s version %s" % (versName, version),
		epilog="""
example: %(prog)s TODO
"""
	)
	parser.add_argument('-i', '--inputs', type=str, metavar='file', action='store', nargs='+', required=True,
		help="input filename(s); may contain wildcards to be expanded using the specified sequence"
	)
	parser.add_argument('-g', '--gprobs-format', action='store_true',
		help="input files are in .gprobs format (3 columns per sample)"
	)
	parser.add_argument('-d', '--dose-format', action='store_true',
		help="input files are in .dose format (1 column per sample)"
	)
	parser.add_argument('-b', '--bgl-files', type=str, metavar='file', action='store', nargs='+', required=True,
		help=".bgl filename(s) from which to replace input 'col.#' headers with the original sample IDs; "+
		"may contain wildcards to be expanded using the specified sequence"
	)
	parser.add_argument('-w', '--wildcard', type=str, metavar='char', action='store', default='#',
		help="wildcard character(s) to expand in input and .bgl filename(s) (default: '#')"
	)
	parser.add_argument('-s', '--sequence', type=str, metavar='list/range', action='store',
		help="sequence with which to expand input and .bgl filename wildcards (default: no wildcard expansion)"
	)
	parser.add_argument('-o', '--output-prefix', type=str, metavar='path', action='store', default='uniq_',
		help="deduplicated output filename prefix (default: 'uniq_')"
	)
	parser.add_argument('-u', '--dupe-prefix', type=str, metavar='path', action='store', default='dupe_',
		help="duplicate comparison output filename prefix (default: 'dupe_')"
	)
	parser.add_argument('-z', '--compressed', type=str, metavar='ext', nargs='?', action='store', default=None, const='',
		help="read compressed input data; with an optional extension, look for input filename(s) plus extension "+
		"from which to read compressed data if available, otherwise use base filename(s) as normal"
	)
	parser.add_argument('-Z', '--output-compressed', action='store_true',
		help="write compressed output data"
	)
	parser.add_argument('-f', '--force', action='store_true',
		help="force overwrite of existing output files"
	)
	parser.add_argument('-v', '--verbose', action='store_true',
		help="print verbose progress messages"
	)
	parser.add_argument('--version', action='version', version=version)
	
	# parse arguments
	args = parser.parse_args()
	
	# validate format flag
	if not (args.gprobs_format or args.dose_format):
		exit("ERROR: one of -g/--gprobs-format or -d/--dose-format must be specified")
	if (args.gprobs_format and args.dose_format):
		exit("ERROR: only one of -g/--gprobs-format or -d/--dose-format may be specified")
	
	# parse and apply wildcard sequence, if specified
	if args.sequence:
		sequence = list()
		for substitute in args.sequence.split(','):
			# try interpreting it as a range
			dash = substitute.find('-',1)
			colon = substitute.find(':',dash+1)
			colon = colon if colon > dash else None
			try:
				start = int(substitute[0:dash])
				stop = int(substitute[dash+1:colon])+1
				step = int(substitute[colon+1:]) if colon else 1
				sequence.extend(str(n) for n in xrange(start,stop,step))
			except ValueError:
				# if int() complains then it's not a numeric range; use it as-is
				sequence.append(substitute)
		#foreach substitute
		
		# do wildcard expansions backwards so the expanding list doesn't break the iteration
		for n in xrange(len(args.inputs)-1,-1,-1):
			if args.wildcard in args.inputs[n]:
				expansion = list()
				for substitute in sequence:
					expansion.append(args.inputs[n].replace(args.wildcard, substitute))
				args.inputs[n:n+1] = expansion
		
		for n in xrange(len(args.bgl_files)-1,-1,-1):
			if args.wildcard in args.bgl_files[n]:
				expansion = list()
				for substitute in sequence:
					expansion.append(args.bgl_files[n].replace(args.wildcard, substitute))
				args.bgl_files[n:n+1] = expansion
	#if args.sequence
	
	# validate file counts
	if len(args.bgl_files) != len(args.inputs):
		exit("ERROR: the number of .bgl files (%d) does not match the number of input files (%d)" % (len(args.bgl_files),len(args.inputs)))
	
	# run
	deduplicate(
		args.inputs,
		args.gprobs_format,
		args.dose_format,
		args.bgl_files,
		args.output_prefix,
		args.dupe_prefix,
		args.compressed,
		args.output_compressed,
		args.force,
		args.verbose
	)
#__main__
