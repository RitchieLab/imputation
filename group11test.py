#!/usr/bin/env python

import sys
import os
import zfile

if len(sys.argv) < 2:
	sys.stderr.write("usage: %s <chromosome>\n" % (sys.argv[0],))
	sys.exit(2)

c = int(sys.argv[1])
gList = (1,2,3,4,5,6,7,8,9,11)

sys.stderr.write("opening chr%d files ...\n" % (c,))
gHeadPath = {}
gHeadFile = {}
for g in gList:
	if g == 11:
		gHeadPath[g] = './eMerge-I_missing/Human660_mapping/chr%d/eMerge660_grp%db37_chr%d_mod.bgl' % (c,g,c)
	else:
		gHeadPath[g] = './eMer_phase1/group%d/Human660_mapping/chr%d/eMerge660grp%db37_chr%d_mod.bgl' % (g,c,g,c)
	gHeadFile[g] = open(gHeadPath[g],'rU')
sys.stderr.write("... OK\n")

# join headers
sys.stderr.write("testing chr%d headers ...\n" % (c,))
sampleFirst = {}
sampleDupes = []
for g in gList:
	# read the header line from each group's file
	words = gHeadFile[g].next().strip().split()
	# append column headers to merged output
	w = 2
	d = 0
	while w < len(words):
		s = int((w / 2) - 1)
		if words[w] in sampleFirst:
			d += 1
			sampleDupes.append( sampleFirst[words[w]]+(g,s) )
		else:
			sampleFirst[words[w]] = (g,s)
		w += 2
	sys.stderr.write("... group %d: %d individuals, %d duplicate samples\n" % (g,((len(words) - 2) / 2),d))
sys.stderr.write("... OK\n")

#output mappings
for dupe in sampleDupes:
	sys.stdout.write("group %d sample %d  ->  group %d sample %d\n" % dupe)
