#!/usr/bin/env python

import sys
import os

if len(sys.argv) < 3:
	sys.stderr.write("usage: %s <snpsPerSeg> <bpsOverlap>\n" % os.path.basename(sys.argv[0]))
	sys.stderr.write("  <snpsPerSeg> = number of markers per segment\n")
	sys.stderr.write("  <bpsOverlap> = number of basepairs of overlap between segments\n")
	exit(1)

snpsPerSeg = int(sys.argv[1])
bpsOverlap = int(sys.argv[2])

sys.stdout.write("chr\tchr.seg\tseg\timp.start.mb\timp.end.mb\tkeep.start.mb\tkeep.end.mb\tsnp.count\tseg.size.mb\n")

for chm in range(1,23):
	snpFile = os.path.join("chr%s" % chm, "chr%s_RHcom.markers" % chm)
	if os.path.exists(snpFile):
		with open(snpFile,'r') as f:
			sys.stderr.write("chr%s ...\n" % chm)
			seg = 1
			num = 0
			p0 = 0
			p1 = 0
			for line in f:
				if num >= snpsPerSeg and p1 >= p0 + 2*bpsOverlap:
					p1 += 1
					o0 = p0 - (bpsOverlap if seg > 1 else 0)
					o1 = p1 + bpsOverlap
					sys.stdout.write("%s\t%s.%s\t%s\t%f\t%f\t%f\t%f\t%d\t%f\n" % (chm,chm,seg,seg,o0/1e6,o1/1e6,p0/1e6,p1/1e6,num,(p1-p0)/1e6))
					seg += 1
					num = 0
					p0 = p1
				num += 1
				p1 = int(line.split()[1])
			if num > 0:
				p1 += 1
				o0 = p0 - (bpsOverlap if seg > 1 else 0)
				o1 = p1
				sys.stdout.write("%s\t%s.%s\t%s\t%f\t%f\t%f\t%f\t%d\t%f\n" % (chm,chm,seg,seg,o0/1e6,o1/1e6,p0/1e6,p1/1e6,num,(p1-p0)/1e6))
			sys.stderr.write("... %d segments\n" % seg)
	else:
		sys.stderr.write("WARNING: %s not found\n" % snpFile)

sys.stderr.write("done.\n")
