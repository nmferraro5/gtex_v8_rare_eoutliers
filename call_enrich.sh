#!/bin/bash

Rscript enrichment/compute_enrichments.R \
        --dir.suffix v8 \
        --outliers.prefix gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.remove \
        --window 10kb_genebody \
	--z.thresh 3 \
	--output.suffix '.v8ciseQTLs.globalOutliers.removed.subset.newParams' \
wait
