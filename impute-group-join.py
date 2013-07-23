#!/usr/bin/env python

import argparse
import collections
import gzip
import os
import sys


if __name__ == "__main__":
	versMaj,versMin,versRev,versDate = 0,7,0,'2013-06-26'
	versStr = "%d.%d.%d (%s)" % (versMaj, versMin, versRev, versDate)
	versDesc = "impute-group-join version %s" % versStr
	
	# define usage
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawDescriptionHelpFormatter,
		description=versDesc,
		epilog="""
example: %(prog)s TODO
"""
	)
	parser.add_argument('-i', '--input', action='append', nargs='+', type=str, metavar='prefix', required=True,
		help="prefix(es) of separate input .dose(.gz) and .gprobs(.gz) files (required)"
	)
	parser.add_argument('-s', '--sample', action='append', nargs='+', type=str, metavar='prefix',
		help="prefix(es) of .bgl(.gz) file(s) for sample header replacement (default: none)"
	)
	parser.add_argument('-f', '--filter', action='store', type=str, metavar='file',
		help="a file listing the sample IDs to retain (default: none)"
	)
	parser.add_argument('-m', '--markers', action='store', type=str, metavar='file',
		help="a file listing the expected order of all markers (default: none)"
	)
	parser.add_argument('-o', '--output', action='store', type=str, metavar='prefix', required=True,
		help="prefix for combined output .dose.gz and .gprobs.gz files"
	)
	parser.add_argument('-d', '--dupes', action='store', type=str, metavar='prefix',
		help="prefix for duplicate sample output .dose.gz and .gprobs.gz files"
	)
	parser.add_argument('--version', action='version', version=versDesc)
	
	# parse arguments
	args = parser.parse_args()
	
	# open input file(s)
	print "finding input .dose(.gz) and .gprobs(.gz) files ..."
	doseFile = list()
	probFile = list()
	for prefixList in args.input:
		for prefix in prefixList:
			if os.path.exists(prefix+'.dose.gz'):
				doseFile.append(gzip.open(prefix+'.dose.gz','rb'))
			elif os.path.exists(prefix+'.dose'):
				doseFile.append(open(prefix+'.dose','rU'))
			else:
				exit("ERROR: %s.dose(.gz) not found" % prefix)
			
			if os.path.exists(prefix+'.gprobs.gz'):
				probFile.append(gzip.open(prefix+'.gprobs.gz','rb'))
			elif os.path.exists(prefix+'.gprobs'):
				probFile.append(open(prefix+'.gprobs','rU'))
			else:
				exit("ERROR: %s.gprobs(.gz) not found" % prefix)
			
			print "  #%d: %s" % (len(doseFile),prefix)
	#foreach args.input
	iRange0 = range(0,len(doseFile))
	iRange1 = range(1,len(doseFile))
	print "... OK: %d sets of input files" % len(doseFile)
	
	# open sample-header file(s) #TODO
	headerFile = None
	if args.sample:
		exit("ERROR: --sample is not yet implemented") #TODO
		headerFile = list()
		print "finding sample header .bgl(.gz) files ..."
		for prefixList in args.sample:
			for prefix in prefixList:
				if os.path.exists(prefix+'.bgl.gz'):
					headerFile.append(gzip.open(prefix+'.bgl.gz','rb'))
				elif os.path.exists(prefix+'.bgl'):
					headerFile.append(open(prefix+'.bgl','rU'))
				else:
					exit("ERROR: %s.bgl(.gz) not found" % prefix)
				
				print "  #%d: %s" % (len(headerFile),prefix)
		#foreach args.sample
		if len(headerFile) != len(doseFile):
			exit("ERROR: input count mismatch (%d .dose/.gprobs files, %d .bgl files)" % (len(doseFile),len(headerFile)))
		print "... OK: %d sample header files" % len(headerFile)
	#if args.sample
	
	# read the marker file, if any
	markerIndex = collections.OrderedDict()
	if args.markers and ((args.markers == '-') or os.path.exists(args.markers)):
		print "reading markers file ..."
		m = 0
		with (sys.stdin if (args.markers == '-') else open(args.markers,'rU')) as markerFile:
			for line in markerFile:
				if not line.startswith('#'):
					marker = line.rstrip("\r\n").split(None,1)[0]
					if marker in markerIndex:
						exit("ERROR: duplicate marker %s" % marker)
					markerIndex[marker] = m
					m += 1
		#with markerFile
		print "... OK: %d markers" % m
	else:
		print "building markers index from input .dose files ..."
		markerList = list()
		markerOrder = dict()
		markerCoverage = dict()
		markerSwap = collections.defaultdict( (lambda: collections.defaultdict(set)) ) # {m1:{m2:{i}}}
		for i in iRange0:
			line = doseFile[i].next()
			if not line.startswith("marker alleleA alleleB"):
				exit("ERROR: invalid header on input .dose file #%d" % (i+1,))
			if i == 0:
				for line in doseFile[i]:
					marker = line[0:line.find(' ')]
					markerOrder[marker] = len(markerList)
					markerList.append(marker)
					markerCoverage[marker] = 0
				#foreach line
				print "  #%d: %d markers" % (i+1,len(markerList))
			else:
				m = mCur = mPrev = mMatch = 0
				for line in doseFile[i]:
					marker = line[0:line.find(' ')]
					try:
						mCur = markerOrder[marker]
						markerCoverage[marker] += 1
						if markerCoverage[marker] > i:
							exit("ERROR: duplicate marker '%s' in input .dose file #%d" % (marker,i+1))
						elif markerCoverage[marker] == i:
							while mCur < mPrev:
								mCur += 1
								if markerCoverage[markerList[mCur]] == i:
									markerSwap[marker][markerList[mCur]].add(i+1)
							mMatch += 1
							mPrev = mCur
						#else:
						#	print "DEBUG: input #%d marker '%s' coverage %d" % (i+1,marker,markerCoverage[marker])
					except KeyError:
						pass
					m += 1
				#foreach line
				print "  #%d: %d markers (%d matching)" % (i+1,m,mMatch)
			#if first input
			doseFile[i].seek(0)
		#foreach input
		i = iRange0[-1]
		
		#print "DEBUG:", markerCoverage
		
		# check for marker swap errors
		iSwaps = None
		for marker1 in markerSwap:
			for marker2 in markerSwap[marker1]:
				if (markerCoverage[marker1] == i) and (markerCoverage[marker2] == i):
					iSwaps = [str(iSwap) for iSwap in sorted(markerSwap[marker1][marker2])]
					print "ERROR: markers '%s' and '%s' order swapped in input .dose file%s #%s" % (marker1,marker2,("" if len(iSwaps) == 1 else "s"),",#".join(iSwaps))
		if iSwaps:
			exit(1)
		
		# compile matched markers
		m = 0
		for marker in markerList:
			if markerCoverage[marker] == i:
				markerIndex[marker] = m
				m += 1
		print "... OK: %d matched markers" % m
		markerList = markerOrder = markerCoverage = markerSwap = None
		
		# write final marker index to file
		if args.markers:
			print "writing markers index file ..."
			with open(args.markers,'wb') as markerFile:
				markerFile.write("\n".join(markerIndex.keys()))
				markerFile.write("\n")
			print "... OK: %d markers written" % len(markerIndex)
	#if args.markers
	
	# initialize buffers
	doseCols = [ None for i in iRange0 ]
	doseUniq = [ list() for i in iRange0 ]
	doseLine = [ None for i in iRange0 ]
	doseSkip = [ 0 for i in iRange0 ]
	probCols = [ None for i in iRange0 ]
	probUniq = [ list() for i in iRange0 ]
	probLine = [ None for i in iRange0 ]
	probSkip = [ 0 for i in iRange0 ]
	
	# read the sample filter file, if any
	sampleFilter = None
	if args.filter:
		sampleFilter = set()
		print "reading sample filter file ..."
		with (sys.stdin if (args.filter == '-') else open(args.filter,'rU')) as filterFile:
			for line in filterFile:
				if not line.startswith('#'):
					sample = line.rstrip("\r\n").split(None,1)[0]
					if sample in sampleFilter:
						exit("ERROR: duplicate sample %s" % sample)
					sampleFilter.add(sample)
		#with markerFile
		print "... OK: %d samples" % (len(sampleFilter),)
	
	# join headers
	print "joining headers ..."
	doseOut = gzip.open(args.output+'.dose.gz', 'wb', compresslevel=6)
	doseOut.write("marker alleleA alleleB")
	probOut = gzip.open(args.output+'.gprobs.gz', 'wb', compresslevel=6)
	probOut.write("marker alleleA alleleB")
	doseDupe = None
	probDupe = None
	sampleFirst = dict()
	sampleDupes = list()
	inputDrop = set()
	for i in iRange0:
		# read sample headers from .bgl headers file if available, otherwise from .dose input file
		nProb = (len(probFile[i].next().rstrip("\r\n").split()) - 3) / 3
		if headerFile: #TODO
			line = headerFile[i].next().rstrip("\r\n").split()
			samples = [sample for s,sample in enumerate(line) if (s > 2) and (s % 2)]
			nHead = len(samples)
			nDose = len(doseFile[i].next().rstrip("\r\n").split()) - 3
			if nProb != nHead or nProb != nDose:
				exit("ERROR: input #%d sample count mismatch (%d in .bgl, %d in .dose, %d in .gprobs)" % (i+1,nHead,nDose,nProb))
		else:
			line = doseFile[i].next().rstrip("\r\n").split()
			samples = line[3:]
			nDose = len(line) - 3
			if nProb != nDose:
				exit("ERROR: input #%d sample count mismatch (%d in .dose, %d in .gprobs)" % (i+1,nDose,nProb))
		
		# identify duplicate samples
		numFilter = numDupe = 0
		for s,sample in enumerate(samples):
			if sampleFilter and (sample not in sampleFilter):
				numFilter += 1
			elif sample in sampleFirst:
				numDupe += 1
				sampleDupes.append(sampleFirst[sample]+(i,s))
				if args.dupes:
					h = " %s(%d/%d)" % (sample,sampleFirst[sample][0],i)
					if not doseDupe:
						doseDupe = gzip.open(args.dupes+'.dose.gz', 'wb', compresslevel=6)
						doseDupe.write("marker alleleA alleleB")
					if not probDupe:
						probDupe = gzip.open(args.dupes+'.gprobs.gz', 'wb', compresslevel=6)
						probDupe.write("marker alleleA alleleB")
					doseDupe.write(h)
					probDupe.write(h*3)
			else:
				sampleFirst[sample] = (i,s)
				doseUniq[i].append(3+s)
				probUniq[i].extend(xrange(3+s*3,3+s*3+3))
				h = " %s" % sample
				doseOut.write(h)
				probOut.write(h*3)
		#foreach samples
		
		# store expected column count
		print "  #%d: %d samples (%d unique, %d duplicate, %d filtered)" % (i+1,len(samples),len(doseUniq[i]),numDupe,numFilter)
		doseCols[i] = 3 + len(samples)
		probCols[i] = 3 + len(samples) * 3
		
		# if any samples will be dropped from this input due to filter or dupe, check if *all* samples were dropped;
		if numFilter or numDupe:
			if not doseUniq[i]:
				doseUniq[i] = False
			if not probUniq[i]:
				probUniq[i] = False
		else:
			doseUniq[i] = True
			probUniq[i] = True
	#foreach input
	doseOut.write("\n")
	probOut.write("\n")
	if doseDupe:
		doseDupe.write("\n")
	if probDupe:
		probDupe.write("\n")
	print "... OK: %d unique samples, %d duplicates" % (len(sampleFirst),len(sampleDupes))
	if sampleDupes and not args.dupes:
		print "WARNING: no duplicate output prefix was specified; duplicate samples will be silently dropped!"
	
	# join lines
	print "joining data ..."
	nextPctP = 10
	nextPctM = int(len(markerIndex) * (nextPctP / 100.0))
	numMatch = 0
	numSkip = 0
	try:
		# initialize input line buffers
		for i in iRange0:
			doseLine[i] = doseFile[i].next().rstrip("\r\n").split()
			probLine[i] = probFile[i].next().rstrip("\r\n").split()
		# join each marker in index order
		for marker,m in markerIndex.iteritems():
			if m > nextPctM:
				print "  ... %d%% ..." % nextPctP
				nextPctP += 10
				nextPctM = int(len(markerIndex) * (nextPctP / 100.0))
			
			# make sure all inputs provide this marker
			match = True
			for i in iRange0:
				while doseLine[i][0] not in markerIndex:
					doseSkip[i] += 1
					doseLine[i] = doseFile[i].next().rstrip("\r\n").split()
				while probLine[i][0] not in markerIndex:
					probSkip[i] += 1
					probLine[i] = probFile[i].next().rstrip("\r\n").split()
				match = match and (doseLine[i][0] == marker) and (probLine[i][0] == marker)
			#foreach input
			if not match:
				numSkip += 1
				# read forward in all files that did contain the matching marker, then move on
				for i in iRange0:
					if doseLine[i][0] == marker:
						doseLine[i] = doseFile[i].next().rstrip("\r\n").split()
					if probLine[i][0] == marker:
						probLine[i] = probFile[i].next().rstrip("\r\n").split()
				#foreach input
				continue
			#if not match
			numMatch += 1
			
			# validate column counts
			for i in iRange0:
				if len(doseLine[i]) != doseCols[i]:
					exit("ERROR: expected %d columns in input .dose file #%d, but found %d for marker '%s'" % (doseCols[i],i+1,len(doseLine[i]),marker))
				if len(probLine[i]) != probCols[i]:
					exit("ERROR: expected %d columns in input .gprobs file #%d, but found %d for marker '%s'" % (probCols[i],i+1,len(probLine[i]),marker))
			#foreach input
			
			# for the first input, store the allele order and then write the data through directly
			a1 = doseLine[0][1]
			a2 = doseLine[0][2]
			if doseUniq[0] == True:
				doseOut.write(" ".join(doseLine[0]))
			elif doseUniq[0] != False:
				doseOut.write("%s %s %s " % (marker,a1,a2))
				doseOut.write(" ".join(doseLine[0][c] for c in doseUniq[0]))
			if probUniq[0] == True:
				probOut.write(" ".join(probLine[0]))
			elif probUniq[0] != False:
				probOut.write("%s %s %s " % (marker,a1,a2))
				probOut.write(" ".join(probLine[0][c] for c in probUniq[0]))
			
			# for other inputs, compare allele order to input 1
			for i in iRange1:
				if doseLine[i][1] == a2 and doseLine[i][2] == a1 and probLine[i][1] == a2 and probLine[i][2] == a1:
					# swap all the gprobs and recalculate doses
					doseLine[i][1] = a1
					doseLine[i][2] = a2
					probLine[i][1] = a1
					probLine[i][2] = a2
					for c in xrange(3,len(probLine[i]),3):
						probLine[i][c],probLine[i][c+2] = probLine[i][c+2],probLine[i][c]
						doseLine[i][2+c/3] = "%f" % (float(probLine[i][c+1])+2*float(probLine[i][c+2]))
					print "  WARNING: swapped allele order for input #%d marker %s" % (i+1,marker)
				elif doseLine[i][1] != a1 or doseLine[i][2] != a2 or probLine[i][1] != a1 or probLine[i][2] != a2:
					exit("ERROR: input #%d marker '%s' allele mismatch (%s/%s expected, %s/%s in .dose, %s/%s in .gprobs)" % (i+1,marker,a1,a2,doseLine[i][1],doseLine[i][2],probLine[i][1],probLine[i][2]))
				if doseUniq[i] == True:
					doseOut.write(" ")
					doseOut.write(" ".join(doseLine[i][3:]))
				elif doseUniq[i] != False:
					doseOut.write(" ")
					doseOut.write(" ".join(doseLine[i][c] for c in doseUniq[i]))
				if probUniq[i] == True:
					probOut.write(" ")
					probOut.write(" ".join(probLine[i][3:]))
				elif probUniq[i] != False:
					probOut.write(" ")
					probOut.write(" ".join(probLine[i][c] for c in probUniq[i]))
			#foreach input
			doseOut.write("\n")
			probOut.write("\n")
			
			# write dupe lines from various inputs, if any
			if sampleDupes and args.dupes:
				doseDupe.write("%s %s %s " % (marker,a1,a2))
				doseDupe.write(" ".join(doseLine[dupe[0]][3+dupe[1]] for dupe in sampleDupes))
				doseDupe.write("\n%s %s %s " % (marker,a1,a2))
				doseDupe.write(" ".join(doseLine[dupe[2]][3+dupe[3]] for dupe in sampleDupes))
				doseDupe.write("\n")
				
				probDupe.write("%s %s %s " % (marker,a1,a2))
				probDupe.write(" ".join(("%s %s %s" % tuple(probLine[dupe[0]][(3+3*dupe[1]):(6+3*dupe[1])])) for dupe in sampleDupes))
				probDupe.write("\n%s %s %s " % (marker,a1,a2))
				probDupe.write(" ".join(("%s %s %s" % tuple(probLine[dupe[2]][(3+3*dupe[3]):(6+3*dupe[3])])) for dupe in sampleDupes))
				probDupe.write("\n")
			#if dupes
			
			# read forward in all files
			for i in iRange0:
				doseLine[i] = doseFile[i].next().rstrip("\r\n").split()
				probLine[i] = probFile[i].next().rstrip("\r\n").split()
			#foreach input
		#foreach marker
	except StopIteration:
		pass
	doseOut.close()
	probOut.close()
	if doseDupe:
		doseDupe.close()
	if probDupe:
		probDupe.close()
	for i in iRange0:
		if doseSkip[i] > 0:
			print "  WARNING: input .dose file #%d had %d extra markers skipped during processing" % (i+1,doseSkip[i])
		n = 0
		try:
			while True:
				doseFile[i].next()
				n += 1
		except StopIteration:
			pass
		doseFile[i].close()
		if n > 0:
			print "  WARNING: input .dose file #%d has %d leftover lines" % (i+1,n)
		
		if probSkip[i] > 0:
			print "  WARNING: input .gprobs file #%d had %d extra markers skipped during processing" % (i+1,probSkip[i])
		n = 0
		try:
			while True:
				probFile[i].next()
				n += 1
		except StopIteration:
			pass
		probFile[i].close()
		if n > 0:
			print "  WARNING: input .gprobs file #%d has %d leftover lines" % (i+1,n)
	#foreach input
	print "... OK: joined %d markers (%d matched, %d incomplete)" % (len(markerIndex),numMatch,numSkip)
#__main__
