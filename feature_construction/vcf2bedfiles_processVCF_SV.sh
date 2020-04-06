#!/bin/bash

set -o nounset

# uses vcftools to extract relevant information from vcf files
# helper script for vcf2bedfiles

usage="usage: vcf2bedfiles_helper_processVCF_SV.sh <SV vcf> <output dir>"
if [ $# -ne 2 ]; then
    echo $usage
    exit
fi

vcfSV=$1
outdir=$2

vcfSVreformat=${outdir}/gtex.lumpy.gs.svscore.high_conf_reformattedMissingGT.vcf
vcfCNV=${outdir}/gtex.lumpy.gs.svscore.high_conf_GSonly.vcf

# individuals to include from the allele frequency calculation
indincl=${RAREDIR}/preprocessing/gtex_2016-01-15_v7_euro_VCFids.txt # only the (123) EA individuals

# fix SV vcf such that the missing genotypes are ./. instead of .
zcat $vcfSV | sed "s%\t\.:%\t./.:%g" | sed 's/""/"/g' > $vcfSVreformat

vcftools --vcf $vcfSVreformat --out ${outdir}/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_GATK_HaplotypeCaller_EAonly_HallLabSV --remove-filtered-all --max-missing-count 10 --keep $indincl --freq &
vcftools --vcf $vcfSVreformat --out ${outdir}/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_GATK_HaplotypeCaller_EAonly_HallLabSV --remove-filtered-all --max-missing-count 10 --keep $indincl --counts &
vcftools --vcf $vcfSVreformat --out ${outdir}/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_GATK_HaplotypeCaller_EAonly_HallLabSV --remove-filtered-all --max-missing-count 10 --keep $indincl --extract-FORMAT-info GT&
vcftools --vcf $vcfSVreformat --out ${outdir}/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_GATK_HaplotypeCaller_EAonly_HallLabSV --remove-filtered-all --max-missing-count 10 --keep $indincl --get-INFO END &

# deal with copy number variants
awk '{if(substr($1,1,1)=="#"){print; next}; if(substr($3,1,2)=="GS"){print}}' $vcfSVreformat > $vcfCNV

vcftools --vcf $vcfCNV --out ${outdir}/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_GATK_HaplotypeCaller_EAonly_HallLabCNV --remove-filtered-all --max-missing-count 10 --keep $indincl --extract-FORMAT-info CN &
vcftools --vcf $vcfCNV --out ${outdir}/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_GATK_HaplotypeCaller_EAonly_HallLabCNV --remove-filtered-all --max-missing-count 10 --keep $indincl --get-INFO END

wait
