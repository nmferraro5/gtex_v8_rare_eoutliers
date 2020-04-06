#!/bin/bash

runEnrich() {
    echo $1
    Rscript enrichment/compute_enrichments_all_types.R --pvalue $1
}

export -f runEnrich

parallel --jobs 4 runEnrich ::: 0.01 0.001 0.0001 0.00001 0.000001 0.0000001 0.00000001 0.000000001
