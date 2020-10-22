#!/bin/bash

### Generate list of top eQTLs for each gene in each tissue

zcat $GTEXv8/eqtl/GTEx_Analysis_v8_eQTL/*.v8.egenes.txt.gz | \
        cut -f14-17 | grep -v variant_pos | sed 's/^X/23/g' | awk 'BEGIN{OFS="\t"}{print $1,$2-1,$2,$3,$4}' | sort -k1,1V | \ #change sorting so it matches vcf i.e. chr1, chr2, etc.
        sed 's/^23/X/g' | uniq > $RAREDIR/preprocessing_v8/gtex_2017-06-05_v8_cis_eQTLs.bed

wait

### Extract these sites from the GTEx v8 VCF using bedtools
zcat $GTEX_WGSv8 | head -4000 | grep '#' > $RAREDIR/preprocessing_v8/gtex_2017-06-05_v8_genotypes_cis_eQTLs.vcf

wait

bedtools intersect -a $RAREDIR/preprocessing_v8/gtex_2017-06-05_v8_cis_eQTLs.bed -b $GTEX_WGSv8 -loj -sorted | \
        grep -v '\-1' | cut -f6- | sort -k1,1 -k2,2n | uniq >> $RAREDIR/preprocessing_v8/gtex_2017-06-05_v8_genotypes_cis_eQTLs.vcf

wait

bgzip $RAREDIR/preprocessing_v8/gtex_2017-06-05_v8_genotypes_cis_eQTLs.vcf

wait

tabix -p vcf $RAREDIR/preprocessing_v8/gtex_2017-06-05_v8_genotypes_cis_eQTLs.vcf.gz

wait
### Convert the cis-eQTL genotypes in VCF format to the number of alternate alleles using VCFTools

vcftools --gzvcf $RAREDIR/preprocessing_v8/gtex_2017-06-05_v8_genotypes_cis_eQTLs.vcf.gz --out $RAREDIR/preprocessing_v8/gtex_2017-06-05_v8_genotypes_cis_eQTLs --012 --maf 0.01

Rscript process_gtex_v8_cis_eqtl_genotypes.R
