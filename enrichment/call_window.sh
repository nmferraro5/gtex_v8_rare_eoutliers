#!/bin/bash

Rscript gene_window_pair_variants.R 100kbp &
wait
Rscript gene_window_pair_variants.R 500kbp &
wait
Rscript gene_window_pair_variants.R 1000kbp &
wait
Rscript gene_window_pair_variants.R 5000kbp &

