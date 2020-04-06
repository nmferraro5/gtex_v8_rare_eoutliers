#!/usr/bin/env python

from __future__ import division
import sys
import argparse
from collections import defaultdict
import pandas as pd
import numpy as np
import os

###
### Command line argument parser
###
parser = argparse.ArgumentParser(description="Build features file with enhancer features.")
parser.add_argument('--features', action='store', default="features.txt", type=str, help='the file to write feature rows to')

### function to remove NAs
def rmNA(values):
        filtered = []
        for val in values:
                if val !="NA":
                        filtered.append(val)
        return filtered

### function to see whether any of the values are non-zero (returns NA if no variants, but returns 0 if not observed)
def get_any(x, colname):
        vals = rmNA(pd.Series(x[colname]))
        if len(vals) == 0: return 'NA'
        return 1 if np.any(map(lambda a: int(a) > 0, vals)) else 0

### function to return chromHMM value for that variant if present
def get_any_ch(x, colname):
        vals = rmNA(pd.Series(x[colname]))
        if len(vals) == 0: return 'NA'
        intvals = [int(a) for a in vals]
        return np.max(intvals) if np.any(map(lambda a: int(a) > 0, vals)) else 0

### function that grabs gene id
def gene_id(x):
	return x['geneID'][0]

###
### Configure input
###
tissues = ['Adipose_Nuclei','Aorta','Colon_Transverse','Sigmoid_Colon','Esophagus','Left_Ventricle','Lung','Pancreas','Peripheral_Blood_Mononuclear_Primary_Cells','Right_Atrium','Muscle_Skeletal','Stomach']
#tissues = ['E027','E062','E063','E065','E066','E068','E069','E071','E073','E074','E075','E076','E079','E080','E086','E095','E096','E097', 'E098','E104','E106','E107','E108','E109','E110','E111','E113','E116','E126']
#tissues = ['Breast','Whole_Blood','Adipose','Aorta','Liver','Esophagus','Lung','Ovary','Pancreas','Intestine','Spleen','LCL','Skin','Brain','Colon','Heart','Muscle_Skeletal','Stomach'] # new enhancers
#columnNames = ['chromosome','TSS0','TSS1','geneID','variantChromosome','Pos0','Pos1','MAF','genotype'] + tissues + ['EAS_AF','AMR_AF','AFR_AF','EUR_AF','SAS_AF']
columnNames = ['chromosome','TSS0','TSS1','geneID','variantChromosome','Pos0','Pos1','MAF','genotype'] + tissues + ['GNOMAD_AF', 'Extra1'] # For HallLabSV
columnNames = ['chromosome','TSS0','TSS1','geneID','variantChromosome','Pos0','Pos1','MAF','genotype'] + tissues + ['GNOMAD_AF']

###
### Main function
###
if __name__ == "__main__":

	### Print command line
	print >> sys.stderr, "[Command]", " ".join(sys.argv)

	### Parse command line arguments
	args = parser.parse_args()

        ### open output file and add header
        featuresFile = open(args.features, 'w')
        featuresFile.write("\t".join(['gene_id'] + tissues) + "\n")

	### Store all the lines in appropriate subgroups per gene
	variants_per_gene = defaultdict(list)

	### Process stdin
	for line in sys.stdin:

		### Breakdown line
		fields = line.strip().split('\t')
                
                ### Use gene ID as key
                key = fields[3]

                ### Store variant data with the gene ID
                variants_per_gene[key].append(fields)
	
	### Iterate over collection of genes
	for gene, variants in variants_per_gene.iteritems():
		### Create a dataframe that we can use logically
		data = pd.DataFrame(variants, columns = columnNames)

		### Collapse variants to create features
		collapsed_features = [gene_id(data)] + [get_any(data, t) for t in tissues]

		featuresFile.write("\t".join([str(s) for s in collapsed_features]) + '\n')

	### Close output file
	featuresFile.close()
