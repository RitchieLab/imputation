#!/usr/bin/env python

import argparse
import collections
import gzip
import itertools
import string
import struct
import sys
import tempfile
import zlib


class zopen(object):
	
	def __init__(self, fileName, splitChar="\n", chunkSize=16*1024):
		self._filePtr = open(fileName,'rb')
		self._splitChar = splitChar
		self._chunkSize = chunkSize
		self._dc = zlib.decompressobj(zlib.MAX_WBITS | 32) # autodetect gzip or zlib header
		self._text = ""
		self._lines = list()
	#__init__()
	
	
	def __del__(self):
		if self._filePtr:
			self._filePtr.close()
	#__del__()
	
	
	def __enter__(self):
		return self
	#__enter__()
	
	
	def __exit__(self, excType, excVal, excTrace):
		pass
	#__exit__()
	
	
	def __iter__(self):
		return self
	#__iter__()
	
	
	def __next__(self):
		# if lines are still cached from the last read, pop one
		if len(self._lines) > 0:
			return self._lines.pop()
		# if there's data left in the source file, read and decompress another chunk
		if self._dc:
			data = self._dc.unused_data
			if data:
				self._dc = zlib.decompressobj(zlib.MAX_WBITS | 32) # autodetect gzip or zlib header
			else:
				data = self._filePtr.read(self._chunkSize)
			if data:
				self._text += self._dc.decompress(data)
				data = None
			else:
				self._text += self._dc.flush()
				self._dc = None
		# if there's no text left, we're done
		if not self._text:
			raise StopIteration
		# split the text into lines
		self._lines = self._text.split(self._splitChar)
		self._text = ""
		# if there's more than one line, store the last to combine with the next chunk
		# (but if there's only one line, and more to read, then keep reading until we get a linebreak)
		if len(self._lines) > 1:
			self._text = self._lines.pop()
		elif self._dc:
			self._text = self._lines.pop()
			self._chunkSize *= 2
			return self.__next__()
		# reverse the remaining lines into a stack and pop one to return
		self._lines.reverse()
		return self._lines.pop()
	#__next__()
	
	
	def next(self):
		return self.__next__()
	#next()
	
	
	def seek(self, offset, whence = 0):
		if offset != 0:
			raise Exception("zfile.seek() does not support offsets != 0")
		self._filePtr.seek(0, whence)
		self._dc.flush()
		self._text = ""
		self._lines = list()
	#seek()
	
#zopen


if __name__ == "__main__":
	versMaj,versMin,versRev,versDate = 0,7,0,'2013-07-08'
	versStr = "%d.%d.%d (%s)" % (versMaj, versMin, versRev, versDate)
	versDesc = "gprobs-to-bed version %s" % versStr
	
	# define usage
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawDescriptionHelpFormatter,
		description=versDesc,
		epilog="""
example: %(prog)s -f my.fam -m my.markers -g my.gprobs.gz -d my.drop -p output
"""
	)
	parser.add_argument('-f', '--fam', type=str, nargs='+', metavar='file', required=True,
		help="input .fam file(s), plain text (required)"
	)
	parser.add_argument('-m', '--markers', type=str, nargs='+', metavar='file', required=True,
		help="input .markers file(s), plain text (required)"
	)
	parser.add_argument('-g', '--gprobs', type=str, nargs='+', metavar='file', required=True,
		help="input .gprobs file(s), compressed (required)"
	)
	parser.add_argument('-d', '--drop', type=str, nargs='+', metavar='file',
		help="input marker drop list file(s), plain text (optional)"
	)
	parser.add_argument('-c', '--chromosome', type=str, metavar='label', default='0',
		help="chromosome to fill in for .map/.bim files (default: 0)"
	)
	parser.add_argument('-t', '--threshold', type=float, metavar='decimal', default=0.51,
		help="minimum probability to call a genotype, in decimal (default: 0.51)"
	)
	parser.add_argument('-p', '--prefix', type=str, metavar='prefix', required=True,
		help="prefix for output files (required)"
	)
	parser.add_argument('--version', action='version', version=versDesc)
	
	# parse arguments
	args = parser.parse_args()
	
	# store some oft-used arguments to save dictionary lookups
	args_chromosome = args.chromosome
	args_threshold = args.threshold
	
	# (1,0,0) -> (1,1)
	# (0,1,0) -> (1,2)
	# (0,0,1) -> (2,2)
	# (0,0,0) -> (0,0)
	
	# read sample file(s)
	sampleLine = dict()
	for famPath in args.fam:
		print "reading sample file '%s' ..." % famPath
		with open(famPath,'rU') as famFile:
			for line in famFile:
				words = line.split(None,2)
				sampleLine['%s_%s' % (words[0],words[1])] = line
		print "... OK: %d samples" % len(sampleLine)
	#for args.fam
	
	# read markers file(s)
	markerPos = dict()
	for markerPath in args.markers:
		print "reading markers file '%s' ..." % markerPath
		with open(markerPath,'rU') as markerFile:
			for line in markerFile:
				words = line.split(None,2)
				markerPos[words[0]] = words[1]
		print "... OK: %d markers" % len(markerPos)
	#for args.markers
	
	# read drop list(s)
	markerDrop = set()
	numDrop = 0
	for dropPath in (args.drop or list()):
		print "reading drop list file '%s' ..." % dropPath
		with open(dropPath,'rU') as dropFile:
			for line in dropFile:
				words = line.split(None,1)
				if words[0] in markerPos:
					markerDrop.add(words[0])
				numDrop += 1
		print "... OK: %d drops (%d applicable)" % (numDrop,len(markerDrop))
	#for args.drop
	
	# read genotype file(s)
	bedFileKeep = open(args.prefix+'.bed','wb')
	bedFileKeep.write("\x6c\x1b") # binary plink file magic number
	bedFileKeep.write("\x01") # SNP-major order (row per snp, columns per sample)
	bimFileKeep = open(args.prefix+'.bim','wb')
	if markerDrop:
		bedFileDrop = open(args.prefix+'.drop.bed','wb')
		bedFileDrop.write("\x6c\x1b") # binary plink file magic number
		bedFileDrop.write("\x01") # SNP-major order (row per snp, columns per sample)
		bimFileDrop = open(args.prefix+'.drop.bim','wb')
	numSamples = 0
	numMarkers = 0
	for genoPath in args.gprobs:
		print "reading genotype file '%s' ..." % genoPath
		with zopen(genoPath) as genoFile:
			header = genoFile.next()
			if not header.startswith("marker alleleA alleleB"):
				print "  ERROR: unexpected file header: %s..." % header[:30]
				sys.exit(1)
			words = header.split()
			nS = (len(words) - 3) / 3
			if not numSamples:
				numSamples = nS
				sWord = range(3,3+3*numSamples,3)
				padding = [0] * (-numSamples % 4)
				byteOut = range(0,(numSamples+len(padding)),4)
				packformat = "%dB"%len(byteOut)
				with open(args.prefix+'.fam','wb') as famFile:
					for s,w in enumerate(sWord):
						line = sampleLine.get(words[w])
						if not line:
							print "  ERROR: unrecognzied sample #%d '%s'" % (s+1,words[w])
							sys.exit(1)
						famFile.write(line)
			elif numSamples != nS:
				print "  ERROR: expected %d samples but found %d" % (numSamples,nS)
				sys.exit(1)
			for line in genoFile:
				words = line.split()
				if words[0] not in markerPos:
					print "  ERROR: unrecognized marker #%d '%s'" % (numMarkers+1,words[0])
					sys.exit(1)
				if words[0] in markerDrop:
					bedFile,bimFile = bedFileDrop,bimFileDrop
				else:
					bedFile,bimFile = bedFileKeep,bimFileKeep
				a1 = a2 = 0
				bedOut = list()
				for w in sWord:
					if float(words[w]) >= args_threshold:
						a1 += 2
						bedOut.append(0b00)
					elif float(words[w+1]) >= args_threshold:
						a1 += 1
						a2 += 1
						bedOut.append(0b10)
					elif float(words[w+2]) >= args_threshold:
						a2 += 2
						bedOut.append(0b11)
					else:
						bedOut.append(0b01)
				#foreach sample
				if a1 > a2:
					if a2 == 0:
						words[2] = "0"
					bimFile.write("%s\t%s\t0\t%s\t%s\t%s\n" % (args_chromosome,words[0],markerPos[words[0]],words[2],words[1]))
					for b,out in enumerate(bedOut):
						if out == 0b00:
							bedOut[b] = 0b11
						elif out == 0b11:
							bedOut[b] = 0b00
				else:
					if a1 == 0:
						words[1] = "0"
					bimFile.write("%s\t%s\t0\t%s\t%s\t%s\n" % (args_chromosome,words[0],markerPos[words[0]],words[1],words[2]))
				bedOut.extend(padding)
				bedFile.write(struct.pack(packformat, *((bedOut[o] | (bedOut[o+1] << 2) | (bedOut[o+2] << 4) | (bedOut[o+3] << 6)) for o in byteOut)))
				numMarkers += 1
			#foreach line in genoFile
		#with genoFile
		print "... OK: %d markers (%d dropped)" % (numMarkers,len(markerDrop))
	#for args.genotype
	bedFileKeep.close()
	bimFileKeep.close()
	if markerDrop:
		bedFileDrop.close()
		bimFileDrop.close()
#__main__

