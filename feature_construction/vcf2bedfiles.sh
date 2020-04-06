#!/bin/bash

# Goes from GTEx VCF file to a bed file for each individual with that individual's variant sites and their allele frequency

set -o nounset

# get current working directory of the script
scriptdir=`dirname \$(readlink -f "\$0")`

vcf=$GTEX_WGSv8
indincl=${RAREDIR}/preprocessing_v8/gtex_2017-06-05_v8_euro_ids.txt
# get individual IDs from subject files - only european ancestry
# first exclude flagged individuals based on wgs data
cat $GTEX_WGSv8_flagged | tail -n +2 | awk '{split($1, id, "-"); print id[1]"-"id[2]}' | \
    grep -v -f - $GTEX_SUBJECTSv8 | awk -F "\t" '$5==3 {print $1}' > $indincl
beddir=${RAREDIR}/features_v8/variantBeds

# create output directory if it doesn't exist
mkdir -p $beddir

# get prefix of output files
fileprefix=`basename $vcf`
fileprefix=${fileprefix%.vcf.gz}
fileprefix=${fileprefix}"_subset"

prefix=${beddir}/${fileprefix}

## actually run things!
#######################
echo "Making variants bed files for $version ..."

# first get vctools to generate useful information
echo "Processing VCF files with vcftools..."
bash ${scriptdir}/vcf2bedfiles_processVCF.sh $vcf $indincl $prefix
echo "Processing VCF files (SNPs/indels) done."

# process SNPs
echo
echo "Processing SNPs..."
bash ${scriptdir}/vcf2bedfiles_processVCFtoolsOutput.sh SNPs $prefix
sleep 5 # so they don't both create the outdir at the same time
echo

# process indels
echo "Processing indels..."
bash ${scriptdir}/vcf2bedfiles_processVCFtoolsOutput.sh indels $prefix
date

# add CADD scores to SNPs
wait
echo
echo "Adding CADD scores to SNPs..."
bash ${scriptdir}/vcf2bedfiles_compile_CADD_scores.sh ${beddir}
echo "Done compiling CADD scores."
date

echo "Gzipping files..."
parallel --jobs 15 gzip ::: ${beddir}/individuals/*.bed
echo "Done gzipping."
