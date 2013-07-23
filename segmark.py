#!/usr/bin/env python

import sys
import argparse


def segmark(inputFiles, chm, segB, segM, mgnB, mgnM, last, outputFile, bases, quiet):
	with (sys.stdout if ((not outputFile) or (outputFile == '-')) else open(outputFile,'w')) as out:
		msg = (sys.stderr if (out == sys.stdout) else sys.stdout)
		
		# enforce seg >= 2*mgn
		if (not quiet):
			msg.write("processing arguments ...\n")
		segB = max(segB, 2*mgnB)
		lstB = max(int(segB * last), 2*mgnB)
		segM = max(segM, 2*mgnM)
		lstM = max(int(segM * last), 2*mgnM)
		if (not quiet):
			msg.write(
"""... OK:
  chr = %s
  segments >= %d bases / %d markers
  margins >= %d bases / %d markers
  last segment >= %d bases / %d markers
""" % (chm or "(unlabeled)", segB, segM, mgnB, mgnM, lstB, lstM)
			)
		
		# read and sort all marker positions
		markerList = [0]
		for inputFile in inputFiles:
			with (sys.stdin if inputFile == '-' else open(inputFile,'r')) as fileObj:
				if (not quiet):
					msg.write("reading marker file '%s' ...\n" % ('<stdin>' if (inputFile == '-') else inputFile))
				n = 0
				for line in fileObj:
					n += 1
					try:
						markerList.append(int(line.split()[0]))
					except ValueError:
						# assume that a value error on the first line means there was a header line
						if n > 1:
							raise
				if (not quiet):
					msg.write("... OK: %d lines\n" % n)
		if (not quiet):
			msg.write("sorting markers ...\n")
		markerList.sort()
		mX = len(markerList) - 1
		bX = markerList[mX]
		if (not quiet):
			msg.write("... OK: %d markers, %d bases\n" % (mX,bX))
		
		# generate segments
		if (not quiet):
			msg.write("segmenting ...\n")
		if chm:
			out.write("#chr\tchr.seg\tseg")
		else:
			out.write("#seg")
		if bases:
			out.write(
				"\tmgnMinB\tsegMinB\tsegMaxB\tmgnMaxB" +
				"\tmgnMinM\tsegMinM\tsegMaxM\tmgnMaxM" +
				"\tmgnSizeB\tsegSizeB\tmgnSizeB" +
				"\tmgnNumM\tsegNumM\tmgnNumM" +
				"\n"
			)
		else:
			out.write(
				"\tmgnMinMB\tsegMinMB\tsegMaxMB\tmgnMaxMB" +
				"\tmgnMinM\tsegMinM\tsegMaxM\tmgnMaxM" +
				"\tmgnSizeMB\tsegSizeMB\tmgnSizeMB" +
				"\tmgnNumM\tsegNumM\tmgnNumM" +
				"\n"
			)
		s = 1
		m = 1
		b = 1
		while (m <= mX) and (b <= bX):
			# calculate margin startpoint
			mgnMinM = max(1, m-mgnM)
			mgnMinB = max(0, b-mgnB)
			if (markerList[mgnMinM] < mgnMinB):
				mgnMinB = markerList[mgnMinM] - ((markerList[mgnMinM] - markerList[max(1,mgnMinM-1)]) // 2)
			else:
				while (mgnMinM > 1) and (markerList[mgnMinM-1] >= mgnMinB):
					mgnMinM -= 1
			
			# calculate segment endpoint
			segMinM = m
			segMinB = b
			segMaxM = min(mX, m+segM-1)
			segMaxB = min(bX, b+segB-1)
			if (markerList[segMaxM] > segMaxB):
				segMaxB = markerList[segMaxM] + ((markerList[min(mX,segMaxM+1)] - markerList[segMaxM]) // 2)
			else:
				while (segMaxM < mX) and (markerList[segMaxM+1] <= segMaxB):
					segMaxM += 1
			
			# if there isn't room for another segment after this one, expand it to fill
			if ((segMaxM + lstM) > mX) or ((segMaxB + lstB) > bX):
				segMaxM = mgnMaxM = mX
				segMaxB = mgnMaxB = bX
			else:
				# calculate margin endpoint
				mgnMaxM = min(mX, segMaxM+mgnM)
				mgnMaxB = min(bX, segMaxB+mgnB)
				if (markerList[mgnMaxM] > mgnMaxB):
					mgnMaxB = markerList[mgnMaxM] + ((markerList[min(mX,mgnMaxM+1)] - markerList[mgnMaxM]) // 2)
				else:
					while (mgnMaxM < mX) and (markerList[mgnMaxM+1] <= mgnMaxB):
						mgnMaxM += 1
			
			# print the segment
			if chm:
				out.write("%s\t%s.%s\t" % (chm,chm,s))
			if bases:
				out.write(
					(
						"%s" +
						"\t%d\t%d\t%d\t%d" +
						"\t%d\t%d\t%d\t%d" +
						"\t%d\t%d\t%d" +
						"\t%d\t%d\t%d" +
						"\n"
					) % (
						s,
						mgnMinB, segMinB, segMaxB, mgnMaxB,
						mgnMinM, segMinM, segMaxM, mgnMaxM,
						(segMinB-mgnMinB), (segMaxB-segMinB+1), (mgnMaxB-segMaxB),
						(segMinM-mgnMinM), (segMaxM-segMinM+1), (mgnMaxM-segMaxM),
					)
				)
			else:
				out.write(
					(
						"%s" +
						"\t%1.6f\t%1.6f\t%1.6f\t%1.6f" +
						"\t%d\t%d\t%d\t%d" +
						"\t%1.6f\t%1.6f\t%1.6f" +
						"\t%d\t%d\t%d" +
						"\n"
					) % (
						s,
						mgnMinB/1e6, segMinB/1e6, segMaxB/1e6, mgnMaxB/1e6,
						mgnMinM,     segMinM,     segMaxM,     mgnMaxM,
						(segMinB-mgnMinB)/1e6, (segMaxB-segMinB+1)/1e6, (mgnMaxB-segMaxB)/1e6,
						(segMinM-mgnMinM),     (segMaxM-segMinM+1),     (mgnMaxM-segMaxM),
					)
				)
			#if bases
			
			s += 1
			m = segMaxM + 1
			b = segMaxB + 1
		#while segments remain
		
		if (not quiet):
			msg.write("... OK: %d segments\n" % (s-1))
	#with outputFile
#segmark()


def segmark_parseSize(arg):
	arg = arg.upper()
	v = b = m = None
	
	if (arg[-2:-1] == 'M'):
		v = int(arg[:-2]) * 1000000
	elif (arg[-2:-1] == 'K'):
		v = int(arg[:-2]) * 1000
	elif (arg[:-1].isdigit()):
		v = int(arg[:-1])
	else:
		sys.stderr.write("ERROR: invalid size '%s', must be numeric\n" % arg)
		exit(2)
	
	if (arg[-1] == 'B'):
		b = v
	elif (arg[-1] == 'M'):
		m = v
	else:
		sys.stderr.write("ERROR: invalid size '%s', must specify 'b' or 'm' suffix\n" % arg)
		exit(2)
	
	return (b,m)
#segmark_parseSize()


if __name__ == "__main__":
	versMaj,versMin,versRev,versDate = 0,2,0,'2012-12-18'
	version = "%d.%d.%d (%s)" % (versMaj, versMin, versRev, versDate)
	
	# define usage
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawDescriptionHelpFormatter,
		description="segmark version %s" % version,
		epilog="""
example: %(prog)s -i chr5.positions -c 5 -s 50mb -m 1mb

The marker positions must be in the first (or only) column of the input
file(s); any extra columns will be ignored.

The positions will be divided into the smallest allowable segments according to
the specified minimum segment and margin sizes.  All sizes must include a
suffix indicating the unit of measure, either 'b' for basepairs or 'm' for
markers; larger unit sizes may be specified with an additional 'k' or 'm'
before the unit suffix.

Normally the last segment could be up to twice the size of the others, since
splitting it would then make it shorter than the minimum; to avoid this, a size
factor may be specified which allows the last segment to be somewhat shorter
than the others.  However, no segment (including the last) may ever be shorter
than twice the margin size.
"""
	)
	parser.add_argument('-i', '--input', action='append', nargs='+', type=str, metavar='file',
		help="list(s) of positions to segment (default: stdin)"
	)
	parser.add_argument('-c', '--chromosome', action='store', type=str, metavar='chr',
		help="chromosome label (default: not used)"
	)
	parser.add_argument('-s', '--segment', action='append', nargs='*', type=str, metavar='size', required=True,
		help="minimum segment size in 'b'ases or 'm'arkers"
	)
	parser.add_argument('-m', '--margin', action='append', nargs='*', type=str, metavar='size', required=True,
		help="minimum margin size in 'b'ases or 'm'arkers"
	)
	parser.add_argument('-l', '--last', action='store', type=float, metavar='decimal', default=1.0,
		help="minimum size factor for the last segment (default: 1.0)"
	)
	parser.add_argument('-o', '--output', action='store', type=str, metavar='file',
		help="output segment map file to write (default: stdout)"
	)
	parser.add_argument('-b', '--bases', action='store_true',
		help="format segment boundaries in single bases instead of the default megabases (default: no)"
	)
	parser.add_argument('-q', '--quiet', action='store_true',
		help="don't print progress messages"
	)
	parser.add_argument('--version', action='version', version=version)
	
	# parse arguments
	args = parser.parse_args()
	
	# collect input file(s)
	inputFiles = []
	if args.input:
		for inputSet in args.input:
			for inputFile in inputSet:
				inputFiles.append(inputFile)
	if len(inputFiles) == 0:
		inputFiles.append('-')
	
	# parse segment size(s)
	segB = segM = 1
	for sizeSet in args.segment:
		for size in sizeSet:
			b,m = segmark_parseSize(size)
			segB = max(segB, b)
			segM = max(segM, m)
	
	# parse margin size(s)
	mgnB = mgnM = 0
	for sizeSet in args.margin:
		for size in sizeSet:
			b,m = segmark_parseSize(size)
			mgnB = max(mgnB, b)
			mgnM = max(mgnM, m)
	
	# run segmentation
	segmark(inputFiles, args.chromosome, segB, segM, mgnB, mgnM, args.last, args.output, args.bases, args.quiet)
#__main__()
