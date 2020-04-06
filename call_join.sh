#!/bin/bash

Rscript enrichment/join_variant_annotations.R --method medz
wait
Rscript enrichment/join_variant_annotations.R --method splicing
wait
Rscript enrichment/join_variant_annotations.R --method ase
