#!/usr/bin/env python

import sys
import os
import gzip
import zfile

if len(sys.argv) < 2:
	print "usage: %s <chromosome>" % (sys.argv[0],)
	exit(0)

c = int(sys.argv[1])
z = int(sys.argv[2]) if len(sys.argv) > 2 else 6
gList = (99,12)

print "opening chr%d files ..." % (c,)
gHeadPath = {}
gHeadFile = {}
gHeadDupe = {}
gProbPath = {}
gProbFile = {}
gProbCols = {}
gProbLine = {}
gProbSkip = {}
gDosePath = {}
gDoseFile = {}
gDoseCols = {}
gDoseLine = {}
gDoseSkip = {}
for g in gList:
	if g == 99:
		gHeadPath[g] = None
	elif g == 12:
		gHeadPath[g] = './eMer_phase1/group%d_2/Human660_mapping/chr%d/eMerge660_grp%d_2b37_chr%d_mod.bgl' % (g,c,g,c)
		gHeadFile[g] = open(gHeadPath[g],'rU')
	else:
		raise Exception
	gHeadDupe[g] = set()
	
	for ext in ('','.gz'):
		if g == 99:
			gProbPath[g] = './results_eMerge-I/chr%d.gprobs%s' % (c,ext)
		elif g == 12:
			gProbPath[g] = './eMer_phase1/group%d_2/Human660_mapping/results/chr%d_grp%d_2.gprobs%s' % (g,c,g,ext)
		else:
			raise Exception
		if os.path.exists(gProbPath[g]):
			gProbFile[g] = zfile.zopen(gProbPath[g]) if (ext == '.gz') else open(gProbPath[g],'rU')
			break
	if g not in gProbFile:
		exit("ERROR: could not find .gprobs or .gprobs.gz file for group %d" % g)
	gProbSkip[g] = 0
	
	for ext in ('','.gz'):
		if g == 99:
			gDosePath[g] = './results_eMerge-I/chr%d.dose%s' % (c,ext)
		elif g == 12:
			gDosePath[g] = './eMer_phase1/group%d_2/Human660_mapping/results/chr%d_grp%d_2.dose%s' % (g,c,g,ext)
		else:
			raise Exception
		if os.path.exists(gDosePath[g]):
			gDoseFile[g] = zfile.zopen(gDosePath[g]) if (ext == '.gz') else open(gDosePath[g],'rU')
			break
	if g not in gDoseFile:
		exit("ERROR: could not find .dose or .dose.gz file for group %d" % g)
	gDoseSkip[g] = 0
if z:
	gProbOut = gzip.GzipFile('./results12/chr%d.gprobs%s' % (c,'.gz' if z else ''), 'wb', compresslevel=z)
	gDoseOut = gzip.GzipFile('./results12/chr%d.dose%s' % (c,'.gz' if z else ''), 'wb', compresslevel=z)
	gProbDupe = gzip.GzipFile('./results12/dupe_chr%d.gprobs%s' % (c,'.gz' if z else ''), 'wb', compresslevel=z)
	gDoseDupe = gzip.GzipFile('./results12/dupe_chr%d.dose%s' % (c,'.gz' if z else ''), 'wb', compresslevel=z)
else:
	gProbOut = open('./results12/chr%d.gprobs%s' % (c,'.gz' if z else ''), 'wb')
	gDoseOut = open('./results12/chr%d.dose%s' % (c,'.gz' if z else ''), 'wb')
	gProbDupe = open('./results12/dupe_chr%d.gprobs%s' % (c,'.gz' if z else ''), 'wb')
	gDoseDupe = open('./results12/dupe_chr%d.dose%s' % (c,'.gz' if z else ''), 'wb')
print "... OK"

# join headers
print "joining chr%d headers ..." % c
gProbOut.write("marker alleleA alleleB")
gDoseOut.write("marker alleleA alleleB")
gProbDupe.write("marker alleleA alleleB")
gDoseDupe.write("marker alleleA alleleB")
sampleFirst = {}
sampleDupes = []
for g in gList:
	# read the header line from each group's file
	if g == 99:
		probs = gProbFile[g].next().strip().split()
		words = [ probs[p] for p in xrange(len(probs)) if (p % 3) != 0]
		gDoseFile[g].next()
	elif g == 12:
		words = gHeadFile[g].next().strip().split()
		gProbFile[g].next()
		gDoseFile[g].next()
	else:
		raise Exception
	# append column headers to merged output
	w = 2
	while w < len(words):
		s = int((w / 2) - 1)
		if words[w] in sampleFirst:
			sampleDupes.append( sampleFirst[words[w]]+(g,s) )
			gHeadDupe[g].add(s)
			h = " %s(%d/%d)" % (words[w],sampleFirst[words[w]][0],g)
			gProbDupe.write(h*3)
			gDoseDupe.write(h)
		else:
			sampleFirst[words[w]] = (g,s)
			h = " %s" % (words[w])
			gProbOut.write(h*3)
			gDoseOut.write(h)
		w += 2
	# store number of expected columns for each group
	gProbCols[g] = 3 + ((len(words) - 2) / 2) * 3
	gDoseCols[g] = 3 + ((len(words) - 2) / 2)
	print "... group %d: %d individuals, %d duplicate samples" % (g,gDoseCols[g]-3,len(gHeadDupe[g]))
gProbOut.write("\n")
gDoseOut.write("\n")
gProbDupe.write("\n")
gDoseDupe.write("\n")
print "... OK"

# join lines
print "joining chr%d results ..." % c
lines = 1
try:
	while True:
		lines += 1
		# read and process data lines for each group
		for g in gList:
			# read each group's next data lines; for groups >1, skip extra markers vs. group 1
			while True:
				gProbLine[g] = gProbFile[g].next().rstrip("\r\n")
				s1 = gProbLine[g].find(' ')
				if g == 99 or gProbLine[g][0:s1] == marker:
					break
				gProbSkip[g] += 1
			while True:
				gDoseLine[g] = gDoseFile[g].next().rstrip("\r\n")
				s1 = gDoseLine[g].find(' ')
				if g == 99 or gDoseLine[g][0:s1] == marker:
					break
				gDoseSkip[g] += 1
			# split lines and validate column count
			gProbLine[g] = gProbLine[g].split()
			gDoseLine[g] = gDoseLine[g].split()
			if len(gProbLine[g]) != gProbCols[g]:
				exit("ERROR: expected %d gprobs columns for group %d, found %d at line %d" % (gProbCols[g],g,len(gProbLine[g]),lines))
			if len(gDoseLine[g]) != gDoseCols[g]:
				exit("ERROR: expected %d dose columns for group %d, found %d at line %d" % (gDoseCols[g],g,len(gDoseLine[g]),lines))
			# for group 99, identify marker and allele order and then write the data through directly
			if g == 99:
				marker = gProbLine[g][0]
				a1 = gProbLine[g][1]
				a2 = gProbLine[g][2]
				gProbOut.write(" ".join(gProbLine[g][w] for w in xrange(0,len(gProbLine[g])) if int((w-3)/3) not in gHeadDupe[g]))
				gDoseOut.write(" ".join(gDoseLine[g][w] for w in xrange(0,len(gDoseLine[g])) if int((w-3)) not in gHeadDupe[g]))
			else:
				# for other groups, compare allele order to group 1
				if gProbLine[g][1] == a2 and gProbLine[g][2] == a1:
					# if the alleles are reversed, swap all the gprobs and recalculate doses
					gProbLine[g][1] = a1
					gProbLine[g][2] = a2
					gDoseLine[g][1] = a1
					gDoseLine[g][2] = a2
					for w in xrange(3,len(gProbLine[g]),3):
						gProbLine[g][w],gProbLine[g][w+2] = gProbLine[g][w+2],gProbLine[g][w]
						gDoseLine[g][2+w/3] = "%f" % (float(gProbLine[g][w+1])+2*float(gProbLine[g][w+2]))
					print "    WARNING: swapped allele order for group %d marker %s (line %d)" % (g,marker,lines)
				if gProbLine[g][1] != a1 or gProbLine[g][2] != a2 or gDoseLine[g][1] != a1 or gDoseLine[g][2] != a2:
					exit("ERROR: group %d marker %s (line %d) allele mismatch: expected=%s/%s, prob=%s/%s, dose=%s/%s" % (g,marker,lines,a1,a2,gProbLine[g][1],gProbLine[g][2],gDoseLine[g][1],gDoseLine[g][2]))
				gProbOut.write(" ")
				gDoseOut.write(" ")
				gProbOut.write(" ".join(gProbLine[g][w] for w in xrange(3,len(gProbLine[g])) if int((w-3)/3) not in gHeadDupe[g]))
				gDoseOut.write(" ".join(gDoseLine[g][w] for w in xrange(3,len(gDoseLine[g])) if int((w-3)) not in gHeadDupe[g]))
			#if group==1
		#foreach group
		gProbOut.write("\n")
		gDoseOut.write("\n")
		
		# write dupe lines from various groups
		gProbDupe.write("%s %s %s " % (marker,a1,a2))
		gDoseDupe.write("%s %s %s " % (marker,a1,a2))
		gProbDupe.write(" ".join(("%s %s %s" % tuple(gProbLine[dupe[0]][(3+3*dupe[1]):(6+3*dupe[1])])) for dupe in sampleDupes))
		gDoseDupe.write(" ".join(gDoseLine[dupe[0]][3+dupe[1]] for dupe in sampleDupes))
		gProbDupe.write("\n%s %s %s " % (marker,a1,a2))
		gDoseDupe.write("\n%s %s %s " % (marker,a1,a2))
		gProbDupe.write(" ".join(("%s %s %s" % tuple(gProbLine[dupe[2]][(3+3*dupe[3]):(6+3*dupe[3])])) for dupe in sampleDupes))
		gDoseDupe.write(" ".join(gDoseLine[dupe[2]][3+dupe[3]] for dupe in sampleDupes))
		gProbDupe.write("\n")
		gDoseDupe.write("\n")
	#forever!
except StopIteration:
	pass
gProbOut.close()
gDoseOut.close()
for g in gList:
	n = 0
	try:
		while True:
			gProbFile[g].next()
			n += 1
	except StopIteration:
		pass
	if gProbSkip[g] > 0:
		print "    WARNING: group %d gprobs file had %d extra lines skipped during processing" % (g,gProbSkip[g])
	if n > 0:
		print "    WARNING: group %d gprobs file has %d leftover lines" % (g,n)
	
	n = 0
	try:
		while True:
			gDoseFile[g].next()
			n += 1
	except StopIteration:
		pass
	if gDoseSkip[g] > 0:
		print "    WARNING: group %d dose file had %d extra lines skipped during processing" % (g,gDoseSkip[g])
	if n > 0:
		print "    WARNING: group %d dose file has %d leftover lines" % (g,n)
print "... OK: %d lines" % lines
