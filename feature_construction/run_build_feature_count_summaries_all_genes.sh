#!/bin/bash

set -o nounset -o pipefail

# get current working directory of the script
scriptdir=`dirname \$(readlink -f "\$0")`


dir=${RAREDIR}/features_v8

indir=${dir}/bySite
indir_sv=${RAREDIR}/features_v7/variantBeds/individuals
outdir=${dir}/byGene
gene=${RAREDIR}/preprocessing_v8/gencode.v26.GRCh38.genes.bed

# Make output directory
mkdir -p ${outdir}

# need to run the features before the counts otherwise the features script won't know to build its sub directory structure

# including the gene body
# 10 kb
bash ${scriptdir}/build_feature_summaries_all_genes.sh $indir ${outdir}/10kb_genebody 10000 $gene
bash ${scriptdir}/build_count_summaries_all_genes.sh $indir ${outdir}/10kb_genebody 10000 $gene

## 200 kb
bash ${scriptdir}/build_feature_summaries_all_genes.sh $indir ${outdir}/200kb_genebody 200000 $gene
bash ${scriptdir}/build_count_summaries_all_genes.sh $indir ${outdir}/200kb_genebody 200000 $gene
