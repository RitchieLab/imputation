#!/usr/bin/env python

import sys
import os
import argparse
import math
import random


def segpheno(phenoPath, filterPath, numGroups, balanceSet, randomize, statsPath, outputPath, outputExt, quiet):
	# read filter data, if any
	filterSet = None
	if filterPath:
		if not quiet:
			sys.stderr.write("reading filter file '%s' ...\n" % filterPath)
		filterSet = set()
		with open(filterPath,'rU') as filterFile:
			# no header on .fam files! pfft
			## validate header
			#filterHeader = filterFile.next().rstrip()
			#words = filterHeader.split()
			#if (len(words) < 2) or (words[1] != "IID"):
			#	sys.stderr.write("ERROR: invalid filter file, expected IID in column 2:\n")
			#	sys.stderr.write(filterHeader)
			#	exit(1)
			#numCols = len(words)
			
			# read IIDs to filter with
			for line in filterFile:
				line = line.rstrip()
				words = line.split()
				#if len(words) != numCols:
				#	sys.stderr.write("ERROR: expected %d columns, found %d\n" % (numCols,len(words)))
				#	exit(1)
				filterSet.add(words[1])
		if not quiet:
			sys.stderr.write("... OK: %d samples\n" % len(filterSet))
	#if filterPath
	
	# read phenotype data
	if not quiet:
		sys.stderr.write("reading phenotype file '%s' ...\n" % phenoPath)
	with open(phenoPath,'rU') as phenoFile:
		# validate header
		phenoHeader = phenoFile.next().rstrip()
		words = phenoHeader.split()
		if (len(words) < 2) or (words[1] != "IID"):
			sys.stderr.write("ERROR: invalid phenotype file, expected IID in column 2:\n")
			sys.stderr.write(phenoHeader)
			exit(1)
		numCols = len(words)
		
		# identify columns to categorize and balance
		colName = {}
		for col,name in enumerate(words):
			if name in balanceSet:
				colName[col] = name
		if len(colName) != len(balanceSet):
			sys.stderr.write("ERROR: expected %d columns to balance, found %d\n" % (len(balanceSet),len(colName)))
			exit(1)
		if statsPath:
			colVal = { col:set() for col in colName }
		
		# read and categorize samples
		catSet = set()
		catRow = {}
		kept = dropped = 0
		for line in phenoFile:
			line = line.rstrip()
			words = line.split()
			if len(words) != numCols:
				sys.stderr.write("ERROR: expected %d columns, found %d\n" % (numCols,len(words)))
				exit(1)
			if (filterSet) and (words[1] not in filterSet):
				dropped += 1
			else:
				kept += 1
				cat = "\t".join([words[col] for col in colName])
				if cat in catSet:
					catRow[cat].append(line)
				else:
					catSet.add(cat)
					if statsPath:
						for col in colName:
							colVal[col].add(words[col])
					catRow[cat] = [line]
		#foreach line
	#with phenoFile
	if not quiet:
		sys.stderr.write("... OK: %d samples in %d categor%s" % (kept,len(catRow),("y" if len(catRow) == 1 else "ies")))
		if filterSet:
			sys.stderr.write(" (%d filtered out)" % dropped)
		sys.stderr.write("\n")
	
	# randomize samples in each category, if requested and needed
	if (randomize) and (len(catSet) > 1):
		if not quiet:
			sys.stderr.write("randomizing categories ...\n")
		# random.shuffle() is faster, but loses randomness on lists longer than
		# 2080 elements due to the Mersenne Twister's period of 2**19937-1;
		# slower alternative: for i in range(0,len(list)): j=truerand(i,len(list)), swap list[i,j]
		for cat in catSet:
			random.shuffle(catRow[cat])
		if not quiet:
			sys.stderr.write("... OK\n")
	
	# distribute samples
	catSize = { cat:len(catRow[cat]) for cat in catSet }
	catNext = { cat:0 for cat in catSet }
	if statsPath:
		grpColValNum = []
	for g in range(0,numGroups):
		grpPath = outputPath + (".grp%d" % (g+1)) + outputExt
		if (not quiet):
			sys.stderr.write("writing group %d to '%s' ...\n" % (g+1,grpPath))
		if statsPath:
			grpColValNum.append( { col:{ val:0 for val in colVal[col] } for col in colName } )
		grpSize = 0
		with open(grpPath,'w') as grpFile:
			grpFile.write(phenoHeader)
			grpFile.write("\n")
			for cat in catSet:
				if catNext[cat] >= catSize[cat]:
					continue
				elif g+1 == numGroups:
					n = catSize[cat] - catNext[cat]
				elif randomize:
					n = 1.0 * catSize[cat] / numGroups
					n = int((n // 1) + (1 if random.random() < (n % 1) else 0))
				else:
					n = int(math.ceil(1.0 * catSize[cat] / numGroups))
				if statsPath:
					vals = cat.split("\t")
				for r in range(catNext[cat], catNext[cat] + n):
					if statsPath:
						v = 0
						for col in colName:
							grpColValNum[g][col][vals[v]] += 1
							v += 1
					grpSize += 1
					grpFile.write(catRow[cat][r])
					grpFile.write("\n")
				catNext[cat] += n
		#with grpFile
		if not quiet:
			sys.stderr.write("... OK: %d samples\n" % grpSize)
	#foreach group
	
	# make sure we wrote out all samples
	for cat in catSet:
		if catNext[cat] != catSize[cat]:
			sys.stderr.write("WARNING: wrote %d of %d samples in a category\n" % (catNext[cat],catSize[cat]))
	
	# output statistics, if requested
	if statsPath:
		if not quiet:
			sys.stderr.write("writing statistics to '%s' ...\n" % ("<stdout>" if statsPath == '-' else statsPath))
		with (sys.stdout if statsPath == '-' else open(statsPath,'w')) as statsFile:
			h1 = ["#"]
			h2 = ["#group"]
			for col in colName:
				for val in colVal[col]:
					h1.append(colName[col])
					h2.append(val)
			statsFile.write("\t".join(h1) + "\n")
			statsFile.write("\t".join(h2) + "\n")
			for g in range(0,numGroups):
				statsFile.write("grp%d" % (g+1))
				for col in colName:
					for val in colVal[col]:
						statsFile.write("\t%d" % grpColValNum[g][col][val])
				statsFile.write("\n")
			if not quiet:
				sys.stderr.write("... OK\n")
	#if statsPath
#segpheno()


if __name__ == "__main__":
	versMaj,versMin,versRev,versDate = 0,1,0,'2011-12-15'
	version = "%d.%d.%d (%s)" % (versMaj, versMin, versRev, versDate)
	
	# define usage
	parser = argparse.ArgumentParser(
		formatter_class=argparse.RawDescriptionHelpFormatter,
		description="segpheno version %s" % version,
		epilog="""
example: %(prog)s -p phenotypes.txt -f filter.fam -n 10 -b SEX RACE -r -o phen

Samples will be read from the input phenotype file and divided into the
specified number of groups.  If a filter is provided, sample IIDs missing from
the filter file will be dropped.  If balance columns are specified, then the
samples will first be categorized according to those columns, and then each
category distributed equally among the output groups.  If statistics are
requested, then the number of samples per group with each observed value of
each balanced column will be reported.
"""
	)
	parser.add_argument('-p', '--pheno', action='store', type=str, metavar='file', required=True,
		help="input phenotype file"
	)
	parser.add_argument('-f', '--filter', action='store', type=str, metavar='file',
		help="if provided, samples from the phenotype file whose IID are not in the filter file will be dropped"
	)
	parser.add_argument('-n', '--num-groups', action='store', type=int, metavar='number', required=True,
		help="number of output groups"
	)
	parser.add_argument('-b', '--balance', action='append', nargs='+', type=str, metavar='column', default=[[]],
		help="categorize samples according to the named columns and distribute categories equally among the output groups"
	)
	parser.add_argument('-r', '--randomize', action='store_true',
		help="randomize group assignments and sample order"
	)
	parser.add_argument('-s', '--stats', action='store', nargs='?', type=str, metavar='file', default=False,
		help="generate statistics about the representation of each balance-column in each group; if no file is specified, write to stdout"
	)
	parser.add_argument('-o', '--output', action='store', type=str, metavar='filename',
		help="output filename prefix (default: input filename without extensions)"
	)
	parser.add_argument('-q', '--quiet', action='store_true',
		help="don't print progress messages"
	)
	parser.add_argument('--version', action='version', version=version)
	
	# parse arguments
	args = parser.parse_args()
	
	# collect balancing columns
	balanceSet = set()
	for colSet in args.balance:
		for col in colSet:
			balanceSet.add(col)
	
	# set default stats output
	if args.stats == None:
		args.stats = '-'
	
	# set default output
	outputPath,outputExt = os.path.splitext(args.pheno)
	outputPath = args.output or outputPath
	outputExt = outputExt or '.txt'
	
	# run join
	segpheno(args.pheno, args.filter, args.num_groups, balanceSet, args.randomize, args.stats, outputPath, outputExt, args.quiet)
#__main__()
