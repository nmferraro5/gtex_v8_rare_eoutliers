#!/bin/bash

set -o nounset -o pipefail

# get current working directory of the script
scriptdir=`dirname \$(readlink -f "\$0")`

## takes as input the version 

dir=${RAREDIR}/features_v8

#indir=${dir}/bySite/HallLabSV
indir=${dir}/bySite
outdir=${dir}/byGene
#outdir=${dir}/byGeneSVs
gene=${RAREDIR}/preprocessing_v8/gencode.v26.GRCh38.genes.bed

# Make output directory
mkdir -p ${outdir}

#indir=${RAREDIR}/features_v7/variantBeds/individuals/HallLabSV_hg38
#indir=${RAREDIR}/features_v7/bySite/HallLabSV
# including the gene body and 500 kb on either side
bash ${scriptdir}/build_enhancer_summaries_all_genes.sh $indir ${outdir}/500kb_genebody 500000 $gene
