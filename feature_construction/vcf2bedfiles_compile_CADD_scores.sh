#!/bin/bash

set -o nounset -o errexit -o pipefail

# get CADD scores for each of the individual SNP bed files

scriptdir=`dirname \$(readlink -f "\$0")`

if [ $# -ne 1 ]; then
    echo "usage: vcf2bedfiles_compile_CADD_scores.sh <bed file directory>"
    exit
fi

indir=${1}/individuals # compatible with outdir of vcf2bedfiles_helper_processVCFtoolsOutput.sh
outdir=${indir}/withCADD

# make output directory if it doesn't exist
if [ ! -d $outdir ]; then
    mkdir $outdir
fi

export indir

date
echo

parallel --jobs 20 "echo {}; cat {} | ${scriptdir}/vcf2bedfiles_compile_CADD_scores.py > {.}.CADD.bed" ::: ${indir}/*_SNPs.bed

# mv all created files to withCADD directory
mv ${indir}/*CADD.bed $outdir

echo
date

echo done
