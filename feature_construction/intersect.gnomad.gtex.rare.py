#!/usr/bin/env python

## Get the number of variants with MAF < 1% and 0.01% per GTEx individual (of european descent).
## The AFs are from gnomad (using the non-finnish european AF).

from __future__ import print_function
import pybedtools
import os

gnomadDir = '/srv/scratch/restricted/GOATs/features_v8/gnomad/'
gtexIndDir = '/srv/scratch/restricted/GOATs/features_v8/variantBeds/individuals/'

## Get list of gnomad files and sort them by chromosome
gnomadFileList = sorted(os.listdir(gnomadDir), key = lambda filename: filename.split('.')[0])
gnomadFileList = [f for f in gnomadFileList if 'all' in f]
## Get list of GTEx files
gtexFileList = os.listdir(gtexIndDir)

## read in gnomad files, combine them, filter for variants with MAF < 0.01 (and MAF < 0.0001)
## produce bed files (in string form) of variants with the desired MAFs
mafStringAll = ''
mafString01 = ''
mafString001 = ''
mafString0001 = ''
for gnomadFile in gnomadFileList:
    gnomad = open(gnomadDir + gnomadFile, 'r')
    dummy = gnomad.readline() # skip header
    for line in gnomad:
        line = line.strip().split()
        ## make AF a float to operate on it (if multiple alt alleles, sum their frequencies)
        try:
		line[4] = sum(map(float, line[4].split(',')))
	except:
		continue
        if line[4] > 0.5:
            line[4] = 1 - line[4] # turn AF into MAF
        ## add string version of the line in bed format to the relevant mafString(s)
        lineString = '\t'.join(['chr' + line[0], str(int(line[1]) - 1), line[1]]) + '\n'
        mafStringAll += lineString
        if line[4] <= 0.01:
            mafString01 += lineString
	    if line[4] <= 0.001:
		mafString001 += lineString
            	if line[4] <= 0.0001:
                	mafString0001 += lineString
    gnomad.close()

## make each of the bed files in string form into bedTool objects
mafAll = pybedtools.BedTool(mafStringAll, from_string = True)
maf01 = pybedtools.BedTool(mafString01, from_string = True)
maf001 = pybedtools.BedTool(mafString001, from_string = True)
maf0001 = pybedtools.BedTool(mafString0001, from_string = True)

## for each GTEx bed file, identify the number of variants not in gnomad
## add this to the number of variants with the predefined MAF in gnomad
for gtexFile in gtexFileList:
    try:
    	gtexBed = pybedtools.BedTool(gtexIndDir + gtexFile)
    	gtexBed01 = maf01.intersect(gtexBed, sorted = True)
    	gtexBed001 = maf001.intersect(gtexBed01, sorted = True)
    	gtexBed0001 = maf0001.intersect(gtexBed001, sorted = True)
    	gtexBedNone = gtexBed.subtract(mafAll)
    	
    	numNone = len(gtexBedNone)
    	num01 = numNone + len(gtexBed01)
    	num001 = numNone + len(gtexBed001)
    	num0001 = numNone + len(gtexBed0001)

    	print('\t'.join([gtexFile, str(num01), str(num001), str(num0001)]))
    except:
        print('\t'.join([gtexFile, 'NA', 'NA', 'NA']))
        continue

