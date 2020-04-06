#!/bin/bash/env

enhDir=/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2017-06-05_v8/references/nygc_roadmap/gtex_v8_roadmap_core15_127/tissue_level_binary/
enhFile=/users/nferraro/data/goats_data/v8_data/enhancer_tissues.txt
outFile=$RAREDIR/preprocessing_v8/roadmap_enhancers.txt
while read ef; do
  #zcat ${enhDir}/${ef}.txt.gz | awk -v tiss=$ef 'BEGIN{OFS="\t"} {print $1,$2,$13,$tiss}' | cut -f1,2,3 >> $outFile
  zcat ${enhDir}/${ef}.txt.gz | cut -f1,2,13 | awk -v tiss=$ef -v OFS="\t" '{print $0,$tiss}' >> $outFile
done <$enhFile

gzip $outFile
