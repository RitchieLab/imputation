#!/usr/bin/env python

import sys
import os
import zlib

if len(sys.argv) <= 5:
	name = os.path.basename(sys.argv[0])
	sys.stderr.write("""
Usage:
  %s <chr> <partsfile> <markerfile> <outputpath> <filetype>
<chr> = chromosome
<partsfile> = segment partition map file used for segmented imputation
<markerfile> = input marker file used for segment imputation
<outputpath> = output directory containing gprobs and dose files
<filetype> = output file type to combine ('dose' or 'gprobs')

Example:
  %s 3 chrom.parts.txt ./chr3_input/chr3.markers ./chr3_output/ dose
will combine the per-segment '.dose.gz' result files in './chr3_output/'
using the chromosome 3 segment margin and keep ranges in 'chrom.parts.txt'
and the original './chr3_input/chr3.markers' file.
""" % (name,name))
	exit(2)


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


# read arguments
chm = sys.argv[1]
partsPath = sys.argv[2]
markerPath = sys.argv[3]
outputPath = sys.argv[4]
outputType = sys.argv[5]

# read segment data
sys.stderr.write("reading segment partition map '%s' ...\n" % partsPath)
segList = [{}] # placeholder segment #0
with open(partsPath,'r') as partsFile:
	header = partsFile.next()
	for line in partsFile:
		words = line.split("\t")
		if words[0] == chm:
			segList.append({
				'mgnStartB': int(float(words[3]) * 1000000),
				'mgnEndB':   int(float(words[4]) * 1000000) - 1,
				'segStartB': int(float(words[5]) * 1000000),
				'segEndB':   int(float(words[6]) * 1000000) - 1
			})
sys.stderr.write("... OK: %d segments\n" % (len(segList)-1))

# read marker data
sys.stderr.write("reading marker file '%s' ...\n" % markerPath)
markerMap = {}
with open(markerPath,'r') as markerFile:
	for line in markerFile:
		words = line.split()
		markerMap[words[0]] = int(words[1])
sys.stderr.write("... OK: %d markers\n" % len(markerMap))

# for each segment, filter and concatenate per-segment outputs
segX = len(segList) - 1
for seg in range(1,segX):
	segStartB = segList[seg]['segStartB']
	segEndB = segList[seg]['segEndB']
	first = True
	resultPath = os.path.join(outputPath, 'output_chr%sseg%d.chr%s_RH_mod.bgl.%s.gz' % (chm,seg,chm,outputType))
	sys.stderr.write("filtering chr%sseg%d %s file '%s' ...\n" % (chm,seg,outputType,resultPath))
	kept = dropped = 0
	for line in zfile(resultPath):
		if first:
			if seg == 1:
				sys.stdout.write(line)
				sys.stdout.write("\n")
			first = False
			continue
		rs = line.split(" ",1)[0]
		pos = markerMap[rs]
		if (seg == 1 or pos >= segStartB) and (seg == segX or pos <= segEndB):
			kept += 1
			sys.stdout.write(line)
			sys.stdout.write("\n")
		else:
			dropped += 1
	sys.stderr.write("... OK: %d kept, %d dropped\n" % (kept,dropped))
#foreach segment
