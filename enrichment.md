Feature enrichment
====================

All paths for created files are relative to `$RAREDIR` unless otherwise noted.

For GTEx v8 samples
-------------------

### Calculate relative risk of rare variants nearby outliers
Collect features for outliers and controls and assess feature enrichment in outliers vs. controls.
Saves the outliers and enrichments in an RData file (`$RAREDIR/data_v8/outliers/enrichments_10kb_genebody_Z3.RData`).
```
Rscript enrichment/compute_enrichments_MEDZ.R \
	--dir.suffix v8 \
	--outliers.prefix gtexV8.outlier.controls. \
	--window 10kb_genebody
```
Relies on `enrichment/enrichment_functions.R`. Can adjust the variant window to other regions. Generates relative risk of rare variants nearby outliers as compared to controls in the given window, in `data_v8/outliers/`.

### Join outlier files with variant annotations
```
Rscript enrichment/join_variant_annotations.R
```
Currently explicitly declares file paths in the script. Relies on variant annotation file, `gtex_v8_rare_GxI_collapsed_feature.tsv`. Generates intersected outlier and variant annotation files in `data_v8/outliers/`.

### Calculate relative risk of variants nearby outliers per type
```
Rscript enrichment/calc_relative_risk_across_types.R
```
Calculates both relative risk of nearby rare variant per category for outliers and association of rare variant status with continuous outlier statistics. Relies on joined outlier and variant annotation files.

### Calculate proportion of variants of each type are nearby outliers
```
Rscript enrichment/get_absolute_risk.R
```
Currently explicitly declares file paths in the script. Relies on variant annotation file, `gtex_v8_rare_GxIxV_variant_annotations.txt.gz`. Generates intersected outlier and variant annotation files in `data_v8/outliers/`.

### Calculate significance of absolute proportions via permutation of outlier status
```
Rscript enrichment/get_permuted_absolute_risk.R
```
Currently explicitly declares file paths in the script. Relies on variant annotation file, `gtex_v8_rare_GxIxV_variant_annotations.txt.gz`. Generates intersected outlier and variant annotation files in `data_v8/outliers/`.



