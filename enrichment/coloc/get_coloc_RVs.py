#!/bin/bash/env python

import pandas as pd
import numpy as np

data_dir = '/users/nferraro/data/goats_data/coloc/'
outliers = pd.read_csv(data_dir + 'coloc_V7outliers.txt', sep='\t')
outliers['UID'] = outliers['Gene'] + '_' + outliers['Ind']
print('Read in outliers')
#var_dir = '/srv/scratch/restricted/GOATs/features_v8/byGene/10kb_genebody/'
var_dir = '/srv/scratch/restricted/GOATs/features_v7/byGene/10kb_genebody/'
snps = pd.read_csv(var_dir + 'all_rare_variants_SNPs.txt', sep='\t',header=None, names=['ID','Gene','Chr', 'Pos','Freq'])
indels = pd.read_csv(var_dir + 'all_rare_variants_indels.txt', sep='\t',header=None, names=['ID','Gene','Chr', 'Pos','Freq'])
print('Read in variants')
print snps.shape
print indels.shape
snps['Gene'] = snps['Gene'].apply(lambda x: x.split('.')[0])
indels['Gene'] = indels['Gene'].apply(lambda x: x.split('.')[0])
snps['UID'] = snps['Gene'] + '_' + snps['ID']
indels['UID'] = indels['Gene'] + '_' + indels['ID']

print list(snps['UID'])[0]

soverlap = [value for value in list(snps['UID']) if value in list(outliers['UID'])]
ioverlap = [value for value in list(indels['UID']) if value in list(outliers['UID'])]

print len(soverlap) 
print len(ioverlap)
