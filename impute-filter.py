#!/usr/bin/env python

import argparse
import sys


if __name__ == "__main__":
	versMaj,versMin,versRev,versDate = 0,5,0,'2013-04-19'
	versStr = "%d.%d.%d (%s)" % (versMaj, versMin, versRev, versDate)
	versDesc = "impute-filter version %s" % versStr
	
	# define usage
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawDescriptionHelpFormatter,
		description="impute-filter",
		epilog="""
example: zcat full.gprobs.gz | %(prog)s -m my.markers | gzip > filtered.gprobs.gz
"""
	)
	parser.add_argument('-m', '--markers', action='store', type=str, metavar='file', required=True,
		help="a file listing the markers that should be kept (required)"
	)
	parser.add_argument('-n', '--numheaders', action='store', type=int, metavar='number', default=0,
		help="number of header rows to always keep (default: 0)"
	)
	
	# parse arguments
	args = parser.parse_args()
	
	# load markers
	markerSet = set()
	sys.stderr.write("reading markers file '%s' ...\n" % args.markers)
	with open(args.markers,'rb') as markerFile:
		for line in markerFile:
			markerSet.add(line.strip())
	sys.stderr.write("... OK: %d markers\n" % len(markerSet))
	
	# filter stdin
	sys.stderr.write("filtering input ...\n")
	for n in xrange(args.numheaders):
		sys.stdout.write(sys.stdin.next())
	n = k = 0
	for line in sys.stdin:
		n += 1
		if line[0:line.find(' ')] in markerSet:
			k += 1
			sys.stdout.write(line)
		# if we know markers will not be repeated,
		# we can bail out a bit early once we've seen them all
	#	if k >= len(markerSet):
	#		break
	sys.stderr.write("... OK: %d lines (%d kept, %d dropped)\n" % (n,k,n-k))
	
