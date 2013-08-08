#!/usr/bin/env python

import argparse
import collections
import gzip
import itertools
import os
import sys


if __name__ == "__main__":
	versMaj,versMin,versRev,versDate = 0,9,0,'2013-08-08'
	versStr = "%d.%d.%d (%s)" % (versMaj, versMin, versRev, versDate)
	versDesc = "impute2-group-join version %s" % versStr
	
	# define usage
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawDescriptionHelpFormatter,
		description=versDesc,
		epilog="""
example: %(prog)s -i a_chr22 b_chr22 -f my.samples -m chr22.markers -o ab_chr22
"""
	)
	parser.add_argument('-i', '--input', action='append', nargs='+', type=str, metavar='prefix', required=True,
		help="prefix(es) of impute2 .phased.sample and .best_guess_haps_imputation.* files to be joined (required)"
	)
	parser.add_argument('-f', '--filter', action='store', type=str, metavar='file',
		help="a file listing the sample IDs to retain (default: none)"
	)
	parser.add_argument('-m', '--markers', action='store', type=str, metavar='file',
		help="a file listing the expected order of all markers (default: none)"
	)
	parser.add_argument('-o', '--output', action='store', type=str, metavar='prefix', required=True,
		help="prefix for joined output and log files"
	)
	parser.add_argument('-d', '--dupes', action='store', type=str, metavar='prefix',
		help="prefix for duplicate sample output files"
	)
	parser.add_argument('--version', action='version', version=versDesc)
	
	# parse arguments
	args = parser.parse_args()
	
	# open input file(s)
	print "finding input files ..."
	prefixList = itertools.chain(*args.input)
	sampleFile = list()
	genoFile = list()
	infoFile = list()
	for prefix in prefixList:
		if os.path.exists(prefix+'.phased.sample.gz'):
			sampleFile.append(gzip.open(prefix+'.phased.sample.gz','rb'))
		elif os.path.exists(prefix+'.phased.sample'):
			sampleFile.append(open(prefix+'.phased.sample','rU'))
		else:
			exit("ERROR: %s.phased.sample(.gz) not found" % prefix)
		
		if os.path.exists(prefix+'.best_guess_haps_imputation.impute2.gz'):
			genoFile.append(gzip.open(prefix+'.best_guess_haps_imputation.impute2.gz','rb'))
		elif os.path.exists(prefix+'.best_guess_haps_imputation.impute2'):
			genoFile.append(open(prefix+'.best_guess_haps_imputation.impute2','rU'))
		else:
			exit("ERROR: %s.best_guess_haps_imputation.impute2(.gz) not found" % prefix)
		
		if os.path.exists(prefix+'.best_guess_haps_imputation.impute2_info.gz'):
			infoFile.append(gzip.open(prefix+'.best_guess_haps_imputation.impute2_info.gz','rb'))
		elif os.path.exists(prefix+'.best_guess_haps_imputation.impute2_info'):
			infoFile.append(open(prefix+'.best_guess_haps_imputation.impute2_info','rU'))
		else:
			exit("ERROR: %s.best_guess_haps_imputation.impute2_info(.gz) not found" % prefix)
		
		print "  #%d: %s" % (len(sampleFile),prefix)
	#foreach prefixList
	iRange0 = range(0,len(sampleFile))
	iRange1 = range(1,len(sampleFile))
	print "... OK: %d sets of input files" % len(sampleFile)
	
	# read the marker file, if any
	markerIndex = collections.OrderedDict()
	if args.markers and ((args.markers == '-') or os.path.exists(args.markers)):
		print "reading markers file ..."
		with (sys.stdin if (args.markers == '-') else open(args.markers,'rU')) as markerFile:
			for line in markerFile:
				if not line.startswith('#'):
					words = line.rstrip("\r\n").split()
					marker = (words[2].lower(), min(words[3],words[4]).lower(), max(words[3],words[4]).lower())
					if marker in markerIndex:
						exit("ERROR: duplicate marker: %s" % (" ".join(words),))
					markerIndex[marker] = words[1]
		#with markerFile
		print "... OK: %d markers" % (len(markerIndex),)
	else:
		print "building markers index from input files ..."
		markerList = list()
		markerOrder = dict()
		markerCoverage = dict()
		markerGeno = dict()
		markerInfo = dict()
		markerLabels = collections.defaultdict(set)
		markerDupe = collections.defaultdict(set)
		markerSwap = collections.defaultdict( lambda:collections.defaultdict(set) ) # {m1:{m2:{i}}}
		
		for i in iRange0:
			m = mCur = mPrev = mMatch = 0
			header = infoFile[i].next().rstrip("\r\n")
			if header != "snp_id rs_id position exp_freq_a1 info certainty type info_type0 concord_type0 r2_type0":
				exit("ERROR: invalid header on info file #%d: %s" % (i+1,header))
			while True:
				# make sure geno/info files agree
				try:
					geno = genoFile[i].next().rstrip("\r\n").split(None,5)[:-1]
				except StopIteration:
					try:
						line = header
						while line == header:
							line = infoFile[i].next().rstrip("\r\n")
						exit("ERROR: input genotype file #%d ended after %d markers, but info file continues" % (i+1,m))
					except StopIteration:
						break
				try:
					line = header
					while line == header:
						line = infoFile[i].next().rstrip("\r\n")
					info = line.split()
				except StopIteration:
					exit("ERROR: input info file #%d ended after %d markers, but genotype file continues" % (i+1,m))
				m += 1
				if (geno[1] != info[1]) or (geno[2] != info[2]):
					exit("ERROR: marker #%d mismatch in input files #%d: '%s %s' vs '%s %s'" % (m,i+1,geno[1],geno[2],info[1],info[2]))
				marker = (geno[2].lower(), min(geno[3],geno[4]).lower(), max(geno[3],geno[4]).lower())
				
				if i == 0:
					# for the first input, just check for duplicates and store metadata
					if marker in markerOrder:
						markerDupe[marker].add(0)
					else:
						markerOrder[marker] = len(markerList)
						markerCoverage[marker] = 1
						markerGeno[marker] = geno
						markerInfo[marker] = info
						for lbl in geno[1].split(';'):
							markerLabels[marker].add(lbl)
					markerList.append(marker)
				else:
					# for subsequent inputs, verify relative order
					try:
						mCur = markerOrder[marker]
						if markerCoverage[marker] > i:
							markerDupe[marker].add(i)
						elif markerCoverage[marker] == i:
							markerCoverage[marker] += 1
							for lbl in geno[1].split(';'):
								markerLabels[marker].add(lbl)
							while mCur < mPrev:
								mCur += 1
								if markerCoverage[markerList[mCur]] > i:
									markerSwap[marker][markerList[mCur]].add(i+1) # +1 here so we can .join() later
							mMatch += 1
							mPrev = mCur
						#else:
						#	print "DEBUG: input #%d marker '%s' coverage %d" % (i+1,marker,markerCoverage[marker])
					except KeyError:
						pass
				#if i
			#while next()
			genoFile[i].seek(0)
			infoFile[i].close()
			if i == 0:
				print "  #%d: %d markers" % (i+1,len(markerList))
			else:
				print "  #%d: %d markers (%d matching)" % (i+1,m,mMatch)
		#foreach input
		
		#print "DEBUG:", markerCoverage
		
		# apply highest RS# labels to all markers
		for marker,labels in markerLabels.iteritems():
			rses = set(int(l[2:]) for l in labels if l.lower().startswith('rs'))
			markerGeno[marker][1] = ('rs%d' % max(rses)) if rses else max(labels)
		
		# check for marker dupe warnings
		if markerDupe:
			print "WARNING: %d markers are duplicated in one or more .impute2(.gz) files" % (len(markerDupe),)
			if args.dupes:
				print "writing duplicate markers to '%s' ..." % (args.dupes+'.markers',)
				with open(args.dupes+'.markers','wb') as dupesFile:
					for marker in sorted(markerDupe, key=markerOrder.get):
						dupesFile.write("%s %s\n" % (" ".join(markerGeno[marker]), " ".join(prefixList[i] for i in sorted(markerDupes[marker]))))
				print "... OK"
		
		# check for marker swap errors
		iSwaps = None
		for marker1 in markerSwap:
			for marker2 in markerSwap[marker1]:
				if (markerCoverage[marker1] == len(iRange0)) and (markerCoverage[marker2] == len(iRange0)):
					iSwaps = [str(iSwap) for iSwap in sorted(markerSwap[marker1][marker2])]
					print "ERROR: marker positions %s and %s order swapped in .impute2(.gz) file(s) #%s" % (marker1[0],marker2[0],",#".join(iSwaps))
		if iSwaps:
			exit(1)
		
		# compile matched markers
		for m,marker in enumerate(markerList):
			if (markerOrder[marker] == m) and (markerCoverage[marker] == len(iRange0)):
				if marker in markerIndex:
					exit("ERROR: duplicate marker: %s" % (" ".join(marker),))
				markerIndex[marker] = markerGeno[marker][1]
		print "... OK: %d matched markers" % (len(markerList),)
		
		# write final marker index to file
		if args.markers:
			print "writing markers index file ..."
			with open(args.markers,'wb') as markerFile:
				markerFile.write("\n".join( " ".join(markerGeno[marker]) for marker in markerIndex.iterkeys() ))
				markerFile.write("\n")
			print "... OK: %d markers written" % len(markerIndex)
		markerList = markerOrder = markerCoverage = markerGeno = markerInfo = markerLabels = markerDupe = markerSwap = None
	#if args.markers
	
	# read the sample filter file, if any
	sampleFilter = None
	if args.filter:
		sampleFilter = set()
		print "reading sample filter file ..."
		with (sys.stdin if (args.filter == '-') else open(args.filter,'rU')) as filterFile:
			for line in filterFile:
				if not line.startswith('#'):
					sample = tuple(line.rstrip("\r\n").lower().split(None,2))[0:2]
					if sample in sampleFilter:
						exit("ERROR: duplicate sample %s" % sample)
					sampleFilter.add(sample)
		#with filterFile
		print "... OK: %d samples" % (len(sampleFilter),)
	
	# initialize buffers
	sampleOut = open(args.output+'.phased.sample', 'wb')
	sampleDupe = None
	genoOut = gzip.open(args.output+'.impute2.gz', 'wb', compresslevel=6)
	genoLog = open(args.output+'.log', 'wb')
	genoLog.write("#chr\tmarker\tpos\tallele1\tallele2\tstatus\tnote\n")
	genoDupe = None
	genoCols = [ None for i in iRange0 ]
	genoUniq = [ list() for i in iRange0 ]
	genoLine = [ None for i in iRange0 ]
	genoMarker = [ None for i in iRange0 ]
	genoSkip = [ 0 for i in iRange0 ]
	
	# check samples
	print "joining samples ..."
	sampleHeader1 = sampleHeader2 = None
	sampleFirst = dict()
	sampleDupes = list()
	inputDrop = set()
	for i in iRange0:
		line = sampleFile[i].next().rstrip("\r\n")
		if not line.startswith("ID_1 ID_2 missing"):
			exit("ERROR: invalid header in .sample input file #%d: %s" % (i+1,line))
		if not sampleHeader1:
			sampleHeader1 = line
			sampleOut.write("%s\n" % sampleHeader1)
		elif line != sampleHeader1:
			exit("ERROR: mismatched header in .sample input file #%d: %s" % (i+1,line))
		
		line = sampleFile[i].next().rstrip("\r\n")
		if not line.startswith("0 0 0"):
			exit("ERROR: invalid subheader in .sample input file #%d: %s" % (i+1,line))
		if not sampleHeader2:
			sampleHeader2 = line
			sampleOut.write("%s\n" % sampleHeader2)
		elif line != sampleHeader2:
			exit("ERROR: mismatched subheader in .sample input file #%d: %s" % (i+1,line))
		
		samples = list()
		for line in sampleFile[i]:
			samples.append(tuple(line.rstrip("\r\n").split()))
		
		# identify duplicate samples
		numFilter = numDupe = 0
		for s,sample in enumerate(samples):
			sampleID = (sample[0].lower(),sample[1].lower())
			if sampleFilter and (sampleID not in sampleFilter):
				numFilter += 1
			elif sampleID in sampleFirst:
				numDupe += 1
				sampleDupes.append(sampleFirst[sampleID]+(i,s))
				if args.dupes:
					if not sampleDupe:
						sampleDupe = open(args.dupes+'.phased.sample', 'wb')
						sampleDupe.write("%s\n" % sampleHeader1)
						sampleDupe.write("%s\n" % sampleHeader2)
					sampleDupe.write("(%d/%d)%s" % (sampleFirst[sampleID][0],i,(" ".join(sample))))
			else:
				sampleFirst[sampleID] = (i,s)
				genoUniq[i].extend(xrange(5+s*3,5+s*3+3))
				sampleOut.write("%s\n" % (" ".join(sample),))
		#foreach samples
		
		# store expected column count
		print "  #%d: %d samples (%d unique, %d duplicate, %d filtered)" % (i+1,len(samples),len(genoUniq[i])/3,numDupe,numFilter)
		genoCols[i] = 5 + len(samples) * 3
		
		# if any samples will be dropped from this input due to filter or dupe, check if *all* samples were dropped;
		if numFilter or numDupe:
			if not genoUniq[i]:
				genoUniq[i] = False
		else:
			genoUniq[i] = True
	#foreach input
	print "... OK: %d unique samples, %d duplicates" % (len(sampleFirst),len(sampleDupes))
	if sampleDupes and not args.dupes:
		print "WARNING: no duplicate output prefix was specified; duplicate samples will be silently dropped!"
	
	# join lines
	print "joining data ..."
	nextPctP = 10
	nextPctM = int(len(markerIndex) * (nextPctP / 100.0))
	numMatch = 0
	numSkip = 0
	markerSkip = set()
	try:
		# initialize input line buffers
		for i in iRange0:
			genoLine[i] = line = genoFile[i].next().rstrip("\r\n").split()
			genoMarker[i] = (line[2].lower(), min(line[3],line[4]).lower(), max(line[3],line[4]).lower())
		# join each marker in index order
		m = 0
		for marker,label in markerIndex.iteritems():
			if m > nextPctM:
				print "  ... %d%% ..." % nextPctP
				nextPctP += 10
				nextPctM = int(len(markerIndex) * (nextPctP / 100.0))
			m += 1
			
			# make sure all inputs provide this marker
			match = True
			for i in iRange0:
				while genoMarker[i] not in markerIndex:
					if genoMarker[i] not in markerSkip:
						markerSkip.add(genoMarker[i])
						genoLog.write("%s\t%s\t%s\t-\tnot matched\n" % (genoLine[i][0], genoLine[i][1], "\t".join(genoMarker[i])))
					genoSkip[i] += 1
					genoLine[i] = line = genoFile[i].next().rstrip("\r\n").split()
					genoMarker[i] = (line[2].lower(), min(line[3],line[4]).lower(), max(line[3],line[4]).lower())
				match = match and (genoMarker[i] == marker)
			#foreach input
			if not match:
				numSkip += 1
				# read forward in all files that did contain the matching marker, then move on
				for i in iRange0:
					if genoMarker[i] == marker:
						if genoMarker[i] not in markerSkip:
							markerSkip.add(genoMarker[i])
							genoLog.write("%s\t%s\t%s\t-\tnot matched\n" % (genoLine[i][0], genoLine[i][1], "\t".join(genoMarker[i])))
						genoLine[i] = line = genoFile[i].next().rstrip("\r\n").split()
						genoMarker[i] = (line[2].lower(), min(line[3],line[4]).lower(), max(line[3],line[4]).lower())
				#foreach input
				continue
			#if not match
			numMatch += 1
			
			# extract marker details, but use the majority label
			chm = genoLine[0][0]
			pos = genoLine[0][2]
			a1 = genoLine[0][3]
			a2 = genoLine[0][4]
			aliases = set(genoLine[i][1] for i in iRange0 if genoLine[i][1] != label)
			if aliases:
				genoLog.write("%s\t%s\t%s\t%s\t%s\t+\t%s\n" % (chm,label,pos,a1,a2,";".join(sorted(aliases))))
			genoLine[0][1] = label
			
			# validate column counts
			for i in iRange0:
				if len(genoLine[i]) != genoCols[i]:
					exit("ERROR: expected %d columns in input .impute2.gz file #%d, but found %d for marker '%s'" % (genoCols[i],i+1,len(genoLine[i]),label))
			#foreach input
			
			# for the first input, store the allele order and then write the data through directly
			if genoUniq[0] == True:
				genoOut.write(" ".join(genoLine[0]))
			elif genoUniq[0] != False:
				genoOut.write("%s %s %s %s %s " % (chm,label,pos,a1,a2))
				genoOut.write(" ".join(genoLine[0][c] for c in genoUniq[0]))
			
			# for other inputs, compare allele order to input 1
			for i in iRange1:
				if genoLine[i][3] == a2 and genoLine[i][4] == a1:
					# swap all the probabilities
					genoLine[i][3] = a1
					genoLine[i][4] = a2
					for c in xrange(5,len(genoLine[i]),3):
						genoLine[i][c],genoLine[i][c+2] = genoLine[i][c+2],genoLine[i][c]
					print "  WARNING: swapped allele order for .impute2(.gz) #%d marker '%s'" % (i+1,label)
				elif genoLine[i][3] != a1 or genoLine[i][4] != a2:
					exit("ERROR: .impute2(.gz) #%d marker '%s' allele mismatch (%s/%s expected, %s/%s found)" % (i+1,label,a1,a2,genoLine[i][3],genoLine[i][4]))
				if genoUniq[i] == True:
					genoOut.write(" ")
					genoOut.write(" ".join(genoLine[i][5:]))
				elif genoUniq[i] != False:
					genoOut.write(" ")
					genoOut.write(" ".join(genoLine[i][c] for c in genoUniq[i]))
			#foreach input
			genoOut.write("\n")
			
			# write dupe lines from various inputs, if any
			if sampleDupes and args.dupes:
				genoDupe.write("%s %s %s %s %s " % (chm,label,pos,a1,a2))
				genoDupe.write(" ".join(("%s %s %s" % tuple(genoLine[dupe[0]][(5+3*dupe[1]):(8+3*dupe[1])])) for dupe in sampleDupes))
				genoDupe.write("\n%s %s %s %s %s " % (chm,label,pos,a1,a2))
				genoDupe.write(" ".join(("%s %s %s" % tuple(genoLine[dupe[2]][(5+3*dupe[3]):(8+3*dupe[3])])) for dupe in sampleDupes))
				genoDupe.write("\n")
			#if dupes
			
			# read forward in all files
			for i in iRange0:
				genoLine[i] = line = genoFile[i].next().rstrip("\r\n").split()
				genoMarker[i] = (line[2].lower(),min(line[3],line[4]).lower(),max(line[3],line[4]).lower())
			#foreach input
		#foreach marker
	except StopIteration:
		pass
	genoOut.close()
	if genoDupe:
		genoDupe.close()
	for i in iRange0:
		if genoSkip[i] > 0:
			print "  WARNING: input .impute2(.gz) file #%d had %d extra markers skipped during processing" % (i+1,genoSkip[i])
		n = 0
		try:
			while True:
				genoFile[i].next()
				n += 1
		except StopIteration:
			pass
		genoFile[i].close()
		if n > 0:
			print "  WARNING: input .impute2(.gz) file #%d has %d leftover lines" % (i+1,n)
	#foreach input
	print "... OK: joined %d markers (%d matched, %d incomplete)" % (len(markerIndex),numMatch,numSkip)
#__main__
