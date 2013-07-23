#!/usr/bin/env python

import sys
import os
import gzip
import zfile

if len(sys.argv) < 2:
	print "usage: %s <chromosome>" % sys.argv[0]
	exit(0)

c = int(sys.argv[1])
print "joining groups for chr: %d" % c

print "opening input files ..."
gBasePath = [
	'./Geisenger_OMNI/Geisenger_b37/group1',
	'./Geisenger_OMNI/Geisenger_b37/group2',
	'./Mayo/Mayo_b37/group1',
	'./Mayo/Mayo_b37/group2',
	'./Mayo/Mayo_b37/group3',
	'./MtSinai_OMNI/group1',
	'./MtSinai_OMNI/group2',
	'./Vanderbilt/group1',
	'./Vanderbilt/group2',
]
gHeadPath = [
	'./Geisenger_OMNI/Geisenger_b37/group1/omni_grp1_chr%d.bgl' % c, 
	'./Geisenger_OMNI/Geisenger_b37/group2/omni_grp2_chr%d.bgl' % c,
	'./Mayo/Mayo_b37/group1/mayo_b37_grp1_chr%d.bgl' % c,
	'./Mayo/Mayo_b37/group2/mayo_b37_grp2_chr%d.bgl' % c,
	'./Mayo/Mayo_b37/group3/mayo_b37_grp3_chr%d.bgl' % c,
	'./MtSinai_OMNI/group1/MtSinai_AA_grp1_chr%d.bgl' % c,
	'./MtSinai_OMNI/group2/MtSinai_AA_grp2_chr%d.bgl' % c,
	'./Vanderbilt/group1/OMNI_VGER_grp1_chr%d.bgl' % c,
	'./Vanderbilt/group2/OMNI_VGER_grp2_chr%d.bgl' % c,
]
gHeadFile = [ None for g in gHeadPath ]
gHeadDupe = [ set() for g in gHeadPath ]
gProbPath = [
	'./Geisenger_OMNI/Geisenger_b37/group1/results/chr%d_omni_grp1.gprobs' % c,
	'./Geisenger_OMNI/Geisenger_b37/group2/results/chr%d_omni_grp2.gprobs' % c,
	'./Mayo/Mayo_b37/group1/results/chr%d_grp1.gprobs' % c,
	'./Mayo/Mayo_b37/group2/results/chr%d_grp2.gprobs' % c,
	'./Mayo/Mayo_b37/group3/results/chr%d_grp3.gprobs' % c,
	'./MtSinai_OMNI/group1/results/chr%d_grp1.gprobs' % c,
	'./MtSinai_OMNI/group2/results/chr%d_grp2.gprobs' % c,
	'./Vanderbilt/group1/results/chr%d_grp1.gprobs' % c,
	'./Vanderbilt/group2/results/chr%d_grp2.gprobs' % c,
]
gProbFile = [ None for g in gProbPath ]
gProbCols = [ None for g in gProbPath ]
gProbLine = [ None for g in gProbPath ]
gProbSkip = [ 0 for g in gProbPath ]
gDosePath = [
	'./Geisenger_OMNI/Geisenger_b37/group1/results/chr%d_omni_grp1.dose' % c,
	'./Geisenger_OMNI/Geisenger_b37/group2/results/chr%d_omni_grp2.dose' % c,
	'./Mayo/Mayo_b37/group1/results/chr%d_grp1.dose' % c,
	'./Mayo/Mayo_b37/group2/results/chr%d_grp2.dose' % c,
	'./Mayo/Mayo_b37/group3/results/chr%d_grp3.dose' % c,
	'./MtSinai_OMNI/group1/results/chr%d_grp1.dose' % c,
	'./MtSinai_OMNI/group2/results/chr%d_grp2.dose' % c,
	'./Vanderbilt/group1/results/chr%d_grp1.dose' % c,
	'./Vanderbilt/group2/results/chr%d_grp2.dose' % c,
]
gDoseFile = [ None for g in gDosePath ]
gDoseCols = [ None for g in gDosePath ]
gDoseLine = [ None for g in gDosePath ]
gDoseSkip = [ 0 for g in gDosePath ]
assert(len(gBasePath)==len(gHeadPath) and len(gBasePath)==len(gProbPath) and len(gBasePath)==len(gDosePath))
groups=xrange(len(gBasePath))
for g in groups:
	print "  group %d: %s" % (g+1,gBasePath[g])
	
	for ext in ('','.gz'):
		if os.path.exists(gHeadPath[g]+ext):
			gHeadFile[g] = zfile.zopen(gHeadPath[g]+ext) if ext else open(gHeadPath[g],'rU')
			break
	if not gHeadFile[g]:
		exit("ERROR: could not find %s with or without .gz extension" % gHeadPath[g])
	
	for ext in ('','.gz'):
		if os.path.exists(gProbPath[g]+ext):
			gProbFile[g] = zfile.zopen(gProbPath[g]+ext) if ext else open(gProbPath[g],'rU')
			break
	if not gProbFile[g]:
		exit("ERROR: could not find %s with or without .gz extension" % gProbPath[g])
	
	for ext in ('','.gz'):
		if os.path.exists(gDosePath[g]+ext):
			gDoseFile[g] = zfile.zopen(gDosePath[g]+ext) if ext else open(gDosePath[g],'rU')
			break
	if not gDoseFile[g]:
		exit("ERROR: could not find %s with or without .gz extension" % gDosePath[g])
print "... OK"

print "opening output files (in ./results/) ..."
gProbOut = gzip.GzipFile('./results/chr%d.gprobs.gz' % c, 'wb', compresslevel=6)
gDoseOut = gzip.GzipFile('./results/chr%d.dose.gz' % c, 'wb', compresslevel=6)
gProbDupe = gzip.GzipFile('./results/dupe_chr%d.gprobs.gz' % c, 'wb', compresslevel=6)
gDoseDupe = gzip.GzipFile('./results/dupe_chr%d.dose.gz' % c, 'wb', compresslevel=6)
print "... OK"

# join headers
print "joining headers ..."
gProbOut.write("marker alleleA alleleB")
gDoseOut.write("marker alleleA alleleB")
gProbDupe.write("marker alleleA alleleB")
gDoseDupe.write("marker alleleA alleleB")
sampleFirst = dict()
sampleDupes = list()
for g in groups:
	# read the header line from each group's file
	words = gHeadFile[g].next().strip().split()
	gProbFile[g].next()
	gDoseFile[g].next()
	# append column headers to merged output
	w = 2
	while w < len(words):
		s = int((w / 2) - 1)
		if words[w] in sampleFirst:
			sampleDupes.append(sampleFirst[words[w]]+(g,s))
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
	print "  group %d: %d individuals, %d duplicates" % (g+1,gDoseCols[g]-3,len(gHeadDupe[g]))
gProbOut.write("\n")
gDoseOut.write("\n")
gProbDupe.write("\n")
gDoseDupe.write("\n")
print "... OK"

# join lines
print "joining dose and gprobs data ..."
lines = 1
try:
	while True:
		lines += 1
		# read and process data lines for each group
		for g in groups:
			# read each group's next data lines; for groups 2-9, skip extra markers vs. group 1
			while True:
				gProbLine[g] = gProbFile[g].next().rstrip("\r\n")
				s1 = gProbLine[g].find(' ')
				if g == groups[0] or gProbLine[g][0:s1] == marker:
					break
				gProbSkip[g] += 1
			while True:
				gDoseLine[g] = gDoseFile[g].next().rstrip("\r\n")
				s1 = gDoseLine[g].find(' ')
				if g == groups[0] or gDoseLine[g][0:s1] == marker:
					break
				gDoseSkip[g] += 1
			# split lines and validate column count
			gProbLine[g] = gProbLine[g].split()
			gDoseLine[g] = gDoseLine[g].split()
			if len(gProbLine[g]) != gProbCols[g]:
				exit("ERROR: expected %d gprobs columns for group %d, found %d at line %d" % (gProbCols[g],g+1,len(gProbLine[g]),lines))
			if len(gDoseLine[g]) != gDoseCols[g]:
				exit("ERROR: expected %d dose columns for group %d, found %d at line %d" % (gDoseCols[g],g+1,len(gDoseLine[g]),lines))
			# for the first group, identify marker and allele order and then write the data through directly
			if g == groups[0]:
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
					print "    WARNING: swapped allele order for group %d marker %s (line %d)" % (g+1,marker,lines)
				if gProbLine[g][1] != a1 or gProbLine[g][2] != a2 or gDoseLine[g][1] != a1 or gDoseLine[g][2] != a2:
					exit("ERROR: group %d marker %s (line %d) allele mismatch: expected=%s/%s, prob=%s/%s, dose=%s/%s" % (g+1,marker,lines,a1,a2,gProbLine[g][1],gProbLine[g][2],gDoseLine[g][1],gDoseLine[g][2]))
				gProbOut.write(" ")
				gDoseOut.write(" ")
				gProbOut.write(" ".join(gProbLine[g][w] for w in xrange(3,len(gProbLine[g])) if int((w-3)/3) not in gHeadDupe[g]))
				gDoseOut.write(" ".join(gDoseLine[g][w] for w in xrange(3,len(gDoseLine[g])) if int((w-3)) not in gHeadDupe[g]))
			#if first group
		#foreach group
		gProbOut.write("\n")
		gDoseOut.write("\n")
		
		# write dupe lines from various groups, if any
		if sampleDupes:
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
for g in groups:
	if gProbSkip[g] > 0:
		print "    WARNING: group %d gprobs file had %d extra lines skipped during processing" % (g+1,gProbSkip[g])
	n = 0
	try:
		while True:
			gProbFile[g].next()
			n += 1
	except StopIteration:
		pass
	#gProbFile[g].close()
	if n > 0:
		print "    WARNING: group %d gprobs file has %d leftover lines" % (g+1,n)
	
	if gDoseSkip[g] > 0:
		print "    WARNING: group %d dose file had %d extra lines skipped during processing" % (g+1,gDoseSkip[g])
	n = 0
	try:
		while True:
			gDoseFile[g].next()
			n += 1
	except StopIteration:
		pass
	#gDoseFile[g].close()
	if n > 0:
		print "    WARNING: group %d dose file has %d leftover lines" % (g+1,n)
print "... OK: %d lines" % lines
