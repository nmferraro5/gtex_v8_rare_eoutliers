#!/bin/bash

# manipulate annotation files in various ways to make them easier to use downstream
# takes in GTEx annotation file from the command line
# Not gzipped, unlike previous

set -o nounset -o errexit -o pipefail

if [ $# -ne 1 ]; then
	echo "Usage: process.reference.files.sh <gtex annotation, gzipped>"
	exit
fi

gtex=$1

rootname=`basename $gtex | sed 's/.gtf//'`
gtexprefix=${RAREDIR}/preprocessing_v8/${rootname}

# gene bed file
cat $gtex | awk 'BEGIN{OFS="\t"}{
if (substr($1,1,1)=="#") {next};
if ($3!="gene") {next};
name=substr($10,2,length($10)-3);
print $1,$4,$5,name
}' > ${gtexprefix}.bed

# gene bed file with 10kb added on either side
cat ${gtexprefix}.bed | awk 'BEGIN{OFS="\t"}{
low = $2-10000;
if (low < 0) {low = 0};
print $1,low,$3+10000,$4
}' > ${gtexprefix}_padded10kb.bed

# genetypes
cat $gtex | awk '{
if($3=="gene" && $1 ~ /^[chr0-9]+$/){
  print substr($10,2,length($10)-3)"\t"substr($14,2,length($14)-3)
}}' > ${gtexprefix}_genetypes_autosomal.txt

# gene bed file with PC and lincRNA only from padded 10kb
grep -E 'lincRNA|protein_coding' ${gtexprefix}_genetypes_autosomal.txt > ${gtexprefix}_genetypes_autosomal_PCandlinc_only.txt
## sort files by gene
sort -k 1 ${gtexprefix}_genetypes_autosomal_PCandlinc_only.txt > tmp.txt && mv tmp.txt ${gtexprefix}_genetypes_autosomal_PCandlinc_only.txt
sort -k 4 ${gtexprefix}_padded10kb.bed > tmp.txt && mv tmp.txt ${gtexprefix}_padded10kb.bed
join --nocheck-order -j1 1 -j2 4 ${gtexprefix}_genetypes_autosomal_PCandlinc_only.txt ${gtexprefix}_padded10kb.bed | awk '{print $2 "\t" $3 "\t" $4 "\t" $1}' > ${gtexprefix}_padded10kb_PCandlinc_only.bed
