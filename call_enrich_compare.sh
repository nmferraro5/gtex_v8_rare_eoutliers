#!/bin/bash


Rscript enrichment/compute_enrichments_MEDZ.R \
    --dir.suffix v8 \
    --outliers.file gtexV8.outlier.controls.peerTop25.medz.txt \
    --window 10kb_genebody \
    --z.thresh 3 \
    --output.suffix 'peerTop25'
wait
Rscript enrichment/compute_enrichments_MEDZ.R \
    --dir.suffix v8 \
    --outliers.file gtexV8.outlier.controls.peerTop50.medz.txt \
    --window 10kb_genebody \
    --z.thresh 3 \
    --output.suffix 'peerTop50'
wait
Rscript enrichment/compute_enrichments_MEDZ.R \
    --dir.suffix v8 \
    --outliers.file gtexV8.outlier.controls.knownOnly.medz.txt \
    --window 10kb_genebody \
    --z.thresh 3 \
    --output.suffix 'knownOnly'
wait
Rscript enrichment/compute_enrichments_MEDZ.R \
    --dir.suffix v8 \
    --outliers.file gtexV8.outlier.controls.nocorrection.medz.txt \
    --window 10kb_genebody \
    --z.thresh 3 \
    --output.suffix 'nocorrection'

