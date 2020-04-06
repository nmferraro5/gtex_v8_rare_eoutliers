#/bin/bash

#Run liftover for all dnase files from hg19 to hg38
#Assumes both hg19tohg38 chain file and liftOver executable have been downloaded

CHAIN=$RAREDIR/preprocessing_v8/hg19ToHg38.over.chain.gz
dtype=dyadic
OUTDIR=$ANNODIR/epigenomicsRoadmap/hg38/$dtype/
for FILE in $(ls $ANNODIR/epigenomicsRoadmap/$dtype/*.bed.gz)
do
  base=$(basename $FILE .bed.gz)
  add='.hg38.bed'
  addun='.hg38.unlifted.bed'
  out=${OUTDIR}$base$add
  outun=${OUTDIR}$base$addun
  cp $FILE ${OUTDIR}
  gunzip ${OUTDIR}${base}.bed.gz
  liftOver ${OUTDIR}${base}.bed $CHAIN $out $outun
  wait
  rm $outun
  rm ${OUTDIR}${base}.bed
  mv $out ${OUTDIR}${base}.bed
  gzip ${OUTDIR}${base}.bed
done
