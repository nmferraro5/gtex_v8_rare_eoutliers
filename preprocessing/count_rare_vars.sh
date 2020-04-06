#!/bin/bash

dataDir=/srv/scratch/restricted/GOATs/features_v8/variantBeds/individuals
outDir=/users/nferraro/data/goats_data/v8_rare/
gnomadDir=/srv/scratch/restricted/GOATs/features_v8/gnomad

#for iF in $(ls $dataDir/*_indels.bed.gz)
#do
#	base=$(basename $iF _indels.bed.gz)
#	maf01=$(zcat $iF | awk 'BEGIN{OFS="\t"}{if($4 < 0.01){print 'NA'}}' | wc -l)
#	maf005=$(zcat $iF | awk 'BEGIN{OFS="\t"}{if($4 < 0.005){print 'NA'}}' | wc -l)
#	maf001=$(zcat $iF | awk 'BEGIN{OFS="\t"}{if($4 < 0.001){print 'NA'}}' | wc -l)
#	echo -e $base $maf01 $maf005 $maf001 | sed 's/ /\t/g' >> ${outDir}all_rare_indels_gtex_v8.txt
#done
#
#wait
#
#for iF in $(ls $dataDir/alreadyDone/*_indels.bed.gz)
#do
#	base=$(basename $iF _indels.bed.gz)
#	maf01=$(zcat $iF | awk 'BEGIN{OFS="\t"}{if($4 < 0.01){print 'NA'}}' | wc -l)
#	maf005=$(zcat $iF | awk 'BEGIN{OFS="\t"}{if($4 < 0.005){print 'NA'}}' | wc -l)
#	maf001=$(zcat $iF | awk 'BEGIN{OFS="\t"}{if($4 < 0.001){print 'NA'}}' | wc -l)
#	echo -e $base $maf01 $maf005 $maf001 | sed 's/ /\t/g' >> ${outDir}all_rare_indels_gtex_v8.txt
#done
#
#wait
#
for iF in $(ls $dataDir/*_SNPs.bed.gz)
do
	base=$(basename $iF _SNPs.bed.gz)
	maf01=$(zcat $iF | awk 'BEGIN{OFS="\t"}{if($4 < 0.01){print 'NA'}}' | wc -l)
	maf005=$(zcat $iF | awk 'BEGIN{OFS="\t"}{if($4 < 0.005){print 'NA'}}' | wc -l)
	maf001=$(zcat $iF | awk 'BEGIN{OFS="\t"}{if($4 < 0.001){print 'NA'}}' | wc -l)
	echo -e $base $maf01 $maf005 $maf001 | sed 's/ /\t/g' >> ${outDir}all_rare_SNPs_gtex_v8.txt
done

wait
