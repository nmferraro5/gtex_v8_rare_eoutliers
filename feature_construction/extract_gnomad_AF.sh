#!/bin/bash

set -o nounset -o errexit -o pipefail

## Extract non-finnish european allele frequencies from gnomad

## only looking at SNVs that are seen in GTEx and that pass filter
##TODO:with MAF <= 0.01

dir=${RAREDIR}/features_v8
outdir=${dir}/gnomad
gtexsnv=${dir}/variantBeds/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased_subset_SNPs_noChr.frq
gtexindel=${dir}/variantBeds/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased_subset_indels_noChr.frq
nproc=5

mkdir -p $outdir

sed 's/^chr\|%$//g' $gtexsnv > ${dir}/variantBeds/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased_subset_SNPs_noChr.frq

for i in {1..22}; do echo $i; done | \
    parallel --jobs $nproc "vcftools --gzvcf ${GNOMAD}/gnomad.genomes.r2.0.2.sites.chr{}.liftover.b38.vcf.gz --out ${outdir}/{}.all --positions $gtexsnv --get-INFO AF_NFE"
for i in {1..22}; do echo $i; done | \
    parallel --jobs $nproc "vcftools --gzvcf ${GNOMAD}/gnomad.genomes.r2.0.2.sites.chr{}.liftover.b38.vcf.gz --out ${outdir}/{}.indel.all --positions $gtexindel --get-INFO AF_NFE"

