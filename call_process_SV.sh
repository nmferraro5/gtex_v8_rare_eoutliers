#!/bin/bash

vcffile='gtex.lumpy.gs.melt.high_conf.vcf.gz'
outdir='/srv/scratch/restricted/GOATs/data_v7/SVs_from_HallLab/'
prefix='GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_GATK_HaplotypeCaller_EAonly'

#bash feature_construction/vcf2bedfiles_processVCF_SV.sh $outdir$vcffile $outdir &> sv.log &
#
#echo "Processing SVs..."
#bash feature_construction/vcf2bedfiles_processVCFtoolsOutput_SV.sh HallLabSV $outdir$prefix
#sleep 5 # so they don't both create the outdir at the same time
#echo
#
#
#svdir='/srv/scratch/restricted/GOATs/features_v7/variantBeds/individuals'
#for svfile in $(ls ${svdir}/*HallLabSV*); do
#	ofile=$(basename $svfile .gz)
#	#zcat $svfile | sort -V -k1,1 -k2,2 > ${svdir}/temp.txt
#	#gunzip $svfile
#	bedtools sort -i $svfile > ${svdir}/temp.txt
#	gzip ${svdir}/temp.txt
#	mv ${svdir}/temp.txt.gz ${svdir}/$ofile.gz
#done
#bash feature_construction/run_build_feature_count_summaries_all_genes.sh v7 &> sv_count.log &
#bash feature_construction/run_add_enhancers_variant_beds.sh v7 5 &> sv_enh.log &
bash feature_construction/run_build_enhancer_summaries_all_genes.sh v7 &> snp_indel_prom.log &
