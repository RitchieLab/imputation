#!/usr/bin/env python

import sys
import argparse
import re
import zlib


def zfile(fileName, splitChar="\n", chunkSize=1*1024*1024):
	dc = zlib.decompressobj(32+zlib.MAX_WBITS) # autodetect gzip or zlib header
	with open(fileName,'rb') as filePtr:
		text = ""
		loop = True
		while loop:
			data = filePtr.read(chunkSize)
			if data:
				text += dc.decompress(data)
				data = None
			else:
				text += dc.flush()
				loop = False
			if text:
				lines = text.split(splitChar)
				i,x = 0,len(lines)-1
				text = lines[x]
				while i < x:
					yield lines[i]
					i += 1
				lines = None
		#while data remains
		if text:
			yield text
#zfile()


def bgljoin(segmarkPath, markerPath, chm, joinSegFile, outputFile, missingCode, quiet, gzip):
	with (sys.stdout if ((not outputFile) or (outputFile == '-')) else open(outputFile,'w')) as out:
		msg = (sys.stderr if (out == sys.stdout) else sys.stdout)
		if gzip:
			outZip = zlib.compressobj()
		
		# read segment data
		if (not quiet):
			msg.write("reading segment map '%s' ...\n" % segmarkPath)
		segList = {}
		with open(segmarkPath,'r') as segmarkFile:
			header = segmarkFile.next()
			if not header.startswith(
				"#chr\tchr.seg\tseg" +
				"\tmgnMinMB\tsegMinMB\tsegMaxMB\tmgnMaxMB" +
				"\tmgnMinM\tsegMinM\tsegMaxM\tmgnMaxM" +
				"\tmgnSizeMB\tsegSizeMB\tmgnSizeMB" +
				"\tmgnNumM\tsegNumM\tmgnNumM"
			):
				msg.write("ERROR: unrecognized header in segment map file:\n")
				msg.write(header)
				exit(1)
			for line in segmarkFile:
				words = line.split("\t")
				if words[0] == chm:
					segList[int(words[2])] = {
						'segMinB': int(float(words[4]) * 1000000),
						'segMaxB': int(float(words[5]) * 1000000),
						'mgnNumM': int(words[14]) + int(words[16]),
						'segNumM': int(words[15])
					}
		if (not quiet):
			msg.write("... OK: %d segments\n" % len(segList))
		
		# verify segment count
		#if len(segList) != len(joinSegFile):
		#	msg.write("ERROR: segment map defines %d segments, but %d data files were provided\n" % (len(segList),len(joinSegFile)))
		#	exit(1)
		
		# read marker data
		if (not quiet):
			msg.write("reading marker file '%s' ...\n" % markerPath)
		markerList = []
		markerMap = {}
		markerA1 = {}
		markerA2 = {}
		with open(markerPath,'r') as markerFile:
			for line in markerFile:
				words = line.split()
				if len(words) < 4:
					msg.write("... ERROR: expected 'marker position allele allele', only found %d columns\n" % len(words))
					exit(1)
				markerList.append(words[0])
				markerMap[words[0]] = int(words[1])
				markerA1[words[0]] = words[2]
				markerA2[words[0]] = words[3]
		if (not quiet):
			msg.write("... OK: %d markers\n" % len(markerMap))
		
		# for each data file, filter out markers in the margins
		header = None
		filler = None
		for seg in sorted(segList.keys()):
			segMinB = segList[seg]['segMinB']
			segMaxB = segList[seg]['segMaxB']
			mgnNumM = segList[seg]['mgnNumM']
			segNumM = segList[seg]['segNumM']
			if seg in joinSegFile:
				if not quiet:
					msg.write("filtering segment %d in '%s' ...\n" % (seg,joinSegFile[seg]))
				firstline = True
				kept = dropped = 0
				for line in zfile(joinSegFile[seg]):
					if (not filler) and (missingCode):
						filler = (" %s" % missingCode) * (line.count(" ") - 2)
					if firstline:
						header = header or line
						if seg == 1:
							if gzip:
								out.write(outZip.compress(line))
								out.write(outZip.compress("\n"))
							else:
								out.write(line)
								out.write("\n")
						firstline = False
					else:
						rs = line.split(" ",1)[0]
						pos = markerMap[rs]
						if (pos >= segMinB) and (pos <= segMaxB):
							kept += 1
							if gzip:
								out.write(outZip.compress(line))
								out.write(outZip.compress("\n"))
							else:
								out.write(line)
								out.write("\n")
						else:
							dropped += 1
				if not quiet:
					if (dropped != mgnNumM) or (kept != segNumM):
						msg.write("... WARNING: %d kept (expected %d), %d dropped (expected %d)\n" % (kept,segNumM,dropped,mgnNumM))
					else:
						msg.write("... OK: %d kept, %d dropped\n" % (kept,dropped))
			elif not missingCode:
				msg.write("ERROR: segment %d is missing\n" % (seg))
				exit(1)
			else:
				if not quiet:
					msg.write("filling in for missing segment %d ...\n" % (seg))
				if (not header) or (not filler):
					msg.write("... ERROR: cannot fill in for first segment, no sample count yet!\n")
					exit(1)
				kept = 0
				for rs in markerList:
					pos = markerMap[rs]
					if (pos >= segMinB) and (pos <= segMaxB):
						kept += 1
						if gzip:
							out.write(outZip.compress("%s %s %s" % (rs,markerA1[rs],markerA2[rs])))
							out.write(outZip.compress(filler))
							out.write(outZip.compress("\n"))
						else:
							out.write("%s %s %s" % (rs,markerA1[rs],markerA2[rs]))
							out.write(filler)
							out.write("\n")
				if not quiet:
					if (dropped != mgnNumM) or (kept != segNumM):
						msg.write("... WARNING: %d filled in (expected %d)\n" % (kept,segNumM))
					else:
						msg.write("... OK:\n")
			#if seg file available
		#foreach seg
		
		if gzip:
			out.write(outZip.flush())
	#with outputFile
#bgljoin()


if __name__ == "__main__":
	versMaj,versMin,versRev,versDate = 0,2,0,'2012-05-08'
	version = "%d.%d.%d (%s)" % (versMaj, versMin, versRev, versDate)
	
	# define usage
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawDescriptionHelpFormatter,
		description="bgljoin version %s" % version,
		epilog="""
example: %(prog)s -s chr5.segmark -m chr5.markers -c 5
         -j output_chr5seg*.bgl.gprobs.gz -M NA

Data from all of the specified Beagle output files will be rejoined using the
original markers and segment map files.  The Beagle files may be provided using
wildcards, as in the example, so long as there are exactly as many files as
the chromosome has segments in the segment map file.  All Beagle data file
names must include 'chr#seg#' so that they can be identified.
"""
	)
	parser.add_argument('-s', '--segmark', action='store', type=str, metavar='file', required=True,
		help="segment map file"
	)
	parser.add_argument('-m', '--markers', action='store', type=str, metavar='file', required=True,
		help="complete original marker file"
	)
	parser.add_argument('-c', '--chromosome', action='store', type=str, metavar='chr', required=True,
		help="chromosome label"
	)
	parser.add_argument('-j', '--join', action='append', nargs='+', type=str, metavar='file', required=True,
		help="Beagle output file(s) to join"
	)
	parser.add_argument('-o', '--output', action='store', type=str, metavar='file',
		help="output file to write (default: stdout)"
	)
	parser.add_argument('-M', '--missing-code', action='store', type=str, metavar='code',
		help="code to use to fill in missing segment data (default: error on missing segment)"
	)
	parser.add_argument('-z', '--gzip', action='store_true',
		help="gzip-compress the output"
	)
	parser.add_argument('-q', '--quiet', action='store_true',
		help="don't print progress messages"
	)
	parser.add_argument('--version', action='version', version=version)
	
	# parse arguments
	args = parser.parse_args()
	
	# collect and identify input files
	joinRegex = re.compile(r'chr[^a-z0-9]?([0-9XYMT]+)[^a-z0-9]?seg[^a-z0-9]?([0-9]+)', re.I)
	joinChrSet = set()
	joinSegFile = {}
	for joinSet in args.join:
		for joinFile in joinSet:
			match = joinRegex.search(joinFile)
			if not match:
				sys.stderr.write("ERROR: could not identify chr/seg for '%s'\n" % joinFile)
				exit(2)
			joinChrSet.add(match.group(1))
			joinSegFile[int(match.group(2))] = joinFile
	if len(joinChrSet) != 1:
		sys.stderr.write("ERROR: filenames suggest %d chromosomes\n" % len(joinChrSet))
		exit(2)
	
	# run join
	bgljoin(args.segmark, args.markers, args.chromosome, joinSegFile, args.output, args.missing_code, args.quiet, args.gzip)
#__main__()
