#!/bin/bash

outpre=gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed
outsuf=.allTissues.v8ciseQTLs.globalOutliers.removed.ase

Rscript enrichment/compute_enrichments_dist.R \
        --dir.suffix v8 \
        --outliers.prefix $outpre \
        --window 1bp_200kb_genebody \
	--z.thresh 3 \
	--output.suffix $outsuf
wait
Rscript enrichment/compute_enrichments_dist.R \
        --dir.suffix v8 \
        --outliers.prefix $outpre \
        --window 200kb_400kb_genebody \
	--z.thresh 3 \
	--output.suffix $outsuf
wait
Rscript enrichment/compute_enrichments_dist.R \
        --dir.suffix v8 \
        --outliers.prefix $outpre \
        --window 400kb_600kb_genebody \
	--z.thresh 3 \
	--output.suffix $outsuf
wait
Rscript enrichment/compute_enrichments_dist.R \
        --dir.suffix v8 \
        --outliers.prefix $outpre \
        --window 600kb_800kb_genebody \
	--z.thresh 3 \
	--output.suffix $outsuf
wait
Rscript enrichment/compute_enrichments_dist.R \
        --dir.suffix v8 \
        --outliers.prefix $outpre \
        --window 800kb_1MB_genebody \
	--z.thresh 3 \
	--output.suffix $outsuf
