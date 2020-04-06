#!/bin/bash

# author: Emily Tsang

# script to put together individual enhancer features for the GTEx WGS samples
# restricts to variants with an AF of 0.25 or less

# INPUT
# * input bed file, see format below
# * 1000 genomes allele frequency file
#   (bed file with 5 additional columns with the AF for each super population [EAS, AMR, AFR, EUR, SAS])
# * output directory

# REQUIRES
# * $KG_AF_INDELS (8 cols, with chr)
# * $KG_AF_SNPS (ditto)
# * $STATEDIR/[prom/dyadic/enh]/*.bed.gz (3 cols, with chr)
# * $RAREDIR/features_v7/er.tissue.map.txt (first column has names of tissues groups)

set -o nounset -o errexit -o pipefail

######################################

# check that there are the correct arguments
if [ $# -ne 3 ]; then
    echo "usage: add_enhancers_variant_beds.sh input_file_name 1KG_AF_bed_file outdir"
    exit
fi

f=$1
kgfile=$2
outdir=$3

fname=`basename $f`
ind=${fname%%.*}
echo $ind
#################################
# input bed file has
# 1* chromosome (with "chr")
# 2* Pos0
# 3* Pos1
# 4* MAF
# 5* genotype (0,1,2)
#
# add the following annotations
# 6-17* presence/absence of enhancer in 12 tissue groups
# 18-22* 1KG allele frequency from 5 super populations (5 columns: EAS, AMR, AFR, EUR, SAS)
##################################

# get prefixes of enhancer files
#tissues=`cat ${RAREDIR}/features_v7/er.tissue.map.txt | tail -n +2 | cut -f1`

# put together command for intersect with variant bed files
# remove last two columns if SNVs
#cmd="zcat $f | awk 'BEGIN{OFS=\"\t\"}{if(NF==7){print \$1,\$2,\$3,\$4,\$5} else {print \$0}}' |"
#for t in $tissues; do
#    cmd="${cmd} intersectBed -sorted -wa -c -a stdin -b ${STATEDIR}/enh/${t}.bed.gz |"
#done
#cmd="${cmd} intersectBed -sorted -wa -wb -loj -a stdin -b ${kgfile} | sed 's/\t\t/\t/g' |"
#cmd="${cmd} awk 'BEGIN{roadmapEnd=17}{for(i=1; i<=roadmapEnd; i++){ printf(\"%s\t\",\$i)}; 
#         for(i=roadmapEnd+4; i<NF; i++){ printf(\"%s\t\",\$i)}; print \$NF}' |"
#cmd="${cmd} sed 's/\t\./\tNA/g' | gzip -c > ${outdir}/${ind}_enhancers.bed.gz"
#
#eval "$cmd"
#
#echo "$fname done"
#
#cmd="zcat $f | awk 'BEGIN{OFS=\"\t\"}{if(NF==7){print \$1,\$2,\$3,\$4,\$5} else {print \$0}}' |"
#for t in $tissues; do
#    cmd="${cmd} intersectBed -wa -c -a stdin -b ${STATEDIR}/enh/${t}.bed.gz |" # removed -sorted
#done
#cmd="${cmd} intersectBed -wa -wb -loj -a stdin -b ${kgfile} | sed 's/\t\t/\t/g' |" # removed -sorted
#cmd="${cmd} awk 'BEGIN{roadmapEnd=17}{for(i=1; i<=roadmapEnd; i++){ printf(\"%s\t\",\$i)}; 
#         for(i=roadmapEnd+4; i<NF; i++){ printf(\"%s\t\",\$i)}; print \$NF}' |"
#cmd="${cmd} sed 's/\t\./\tNA/g' | gzip -c > ${outdir}/${ind}_enhancers.bed.gz"
#
#eval "$cmd"


cmd="zcat $f | awk 'BEGIN{OFS=\"\t\"}{if(NF==7){print \$1,\$2,\$3,\$4,\$5} else {print \$0}}' |"
#cmd="${cmd} intersectBed -wa -wb -a stdin -b ${RAREDIR}/preprocessing_v8/gtex_v8.core15.127.filtered.enhancers.binary.txt.gz |" # removed -sorted
cmd="${cmd} intersectBed -wa -wb -a stdin -b ${RAREDIR}/preprocessing_v8/gtex_v8.core15.127.filtered.enhancers.txt.gz |" # removed -sorted
cmd="${cmd} intersectBed -wa -wb -loj -a stdin -b ${kgfile} | sed 's/\t\t/\t/g' |" # removed -sorted
cmd="${cmd} sed 's/\t\./\tNA/g' | cut -f1-4,9-37,41 | gzip -c > ${outdir}/${ind}_chromHMM.bed.gz"
#cmd="${cmd} sed 's/\t\./\tNA/g' | cut -f1-5,9-26,30 | gzip -c > ${outdir}/${ind}_chromHMM.bed.gz" # used in enhancersNew

#cmd="zcat $f | awk 'BEGIN{OFS=\"\t\"}{if(NF==7){print \$1,\$2,\$3,\$4,\$5} else {print \$0}}' |"
#cmd="${cmd} intersectBed -wa -c -a stdin -b ${RAREDIR}/preprocessing_v8/gtex_v8.core15.127.filtered.enhancers.binary.txt |" # removed -sorted
#cmd="${cmd} intersectBed -wa -wb -loj -a stdin -b ${kgfile} | sed 's/\t\t/\t/g' |" # removed -sorted
#cmd="${cmd} awk 'BEGIN{roadmapEnd=17}{for(i=1; i<=roadmapEnd; i++){ printf(\"%s\t\",\$i)}; 
#         for(i=roadmapEnd+4; i<NF; i++){ printf(\"%s\t\",\$i)}; print \$NF}' |"
#cmd="${cmd} sed 's/\t\./\tNA/g' | gzip -c > ${outdir}/${ind}_enhancersNew.bed.gz"

eval "$cmd"
echo "$fname done"
