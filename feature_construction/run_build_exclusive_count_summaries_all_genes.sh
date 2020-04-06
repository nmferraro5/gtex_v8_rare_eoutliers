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

# 0 to 200kb 
bash ${scriptdir}/build_count_exclusive_summaries_all_genes.sh $VERSION $indir ${outdir}/10kb_genebody 0 10000 $gene
wait
bash ${scriptdir}/build_count_exclusive_summaries_all_genes.sh $VERSION $indir ${outdir}/10kb_50kb_genebody 10000 50000 $gene
wait
bash ${scriptdir}/build_count_exclusive_summaries_all_genes.sh $VERSION $indir ${outdir}/50kb_100kb_genebody 50000 100000 $gene
wait
bash ${scriptdir}/build_count_exclusive_summaries_all_genes.sh $VERSION $indir ${outdir}/100kb_200kb_genebody 100000 200000 $gene
wait
## 10kb to 100kb
bash ${scriptdir}/build_count_exclusive_summaries_all_genes.sh $VERSION $indir ${outdir}/200kb_350kb_genebody 200000 350000 $gene
wait
#
## 100kb to 200kb
bash ${scriptdir}/build_count_exclusive_summaries_all_genes.sh $VERSION $indir ${outdir}/350kb_500kb_genebody 350000 500000 $gene
wait
#
## 200kb to 500kb
bash ${scriptdir}/build_count_exclusive_summaries_all_genes.sh $VERSION $indir ${outdir}/200kb_500kb_genebody 200000 500000 $gene
wait

## 500kb to 750kb
bash ${scriptdir}/build_count_exclusive_summaries_all_genes.sh $VERSION $indir ${outdir}/500kb_750kb_genebody 500000 750000 $gene
wait

# 750kb to 1MB
bash ${scriptdir}/build_count_exclusive_summaries_all_genes.sh $VERSION $indir ${outdir}/750kb_1MB_genebody 750000 1000000 $gene
