Data correction, transformation, and organization
=================================================

The list of created files are given as relative paths under `$RAREDIR`.

For GTEx v8 samples
-------------------

### Split data by tissue
```
OUT=${RAREDIR}/preprocessing_v8/PEER_v8
GTEX=${GTEX_RNAv8}/GTEx_Analysis_2017-06-05_v8_RSEMv1.3.0_gene_tpm.gct.gz
SAMPLE=${RAREDIR}/preprocessing_v8/gtex_2017-06-05_v8_samples_tissues.txt
END='.tpm.txt'

python preprocessing/split_expr_by_tissues.py --GTEX $GTEX --OUT $OUT --SAMPLE $SAMPLE --END $END

GTEX=${GTEX_RNAv8}/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz
END='.reads.txt'

python preprocessing/split_expr_by_tissues.py --GTEX $GTEX --OUT $OUT --SAMPLE $SAMPLE --END $END
```
Creates one file with read counts and one with tpm per tissue in `preprocessing_v8` folder

### Transforming data prior to PEER correction
```
Rscript preprocessing/preprocess_expr.R
```
The generated mapping file excludes flagged individuals and samples.
Creates `preprocessing_v8/gtex_2017-06-05_v8_samples_tissues.txt` and some files under `preprocessing_v8/PEER_v8/`.

### Generate list of top eQTLs for each gene in each tissue, extract from VCF, convert to number alternative alleles
```
bash get_eqtl_genotypes.sh
```
Generates several intermediate files in `preprocessing_v8` and relies on `process_gtex_v8_cis_eqtl_genotypes.R` to generate final `gtex_2017-06-05_v8_genotypes_cis_eQTLs_012_processed.txt` in `preprocessing_v8`

### Actually run PEER correction and compute residuals
```
bash preprocessing/correction/calculate_PEER.sh
```
Relies on `preprocessing/correction/calculate_PEER_factors.R` and `preprocessing/correction/calculate_PEER_residuals.R`.

Creates a file for each tissue  under `preprocessing/PEER_v8/` with scaled and corrected log2(tpm) values.

### Combine PEER-corrected data into a single flat file (and compress output file)
```
python gather_filter_normalized_expression.py 
gzip ${RAREDIR}/preprocessing_v8/gtex_2017-06-05_normalized_expression.txt
```

Creates and compresses `preprocessing_v8/gtex_2017-06-05_normalized_expression.txt.gz`.

### Generate files with data on what tissues are available per individual
```
bash get_tissue_by_individual.sh
```
Generates `gtex_2017-06-05_tissues_all_normalized_samples.txt` and `gtex_2017-06-05_tissues_all_normalized_samples.txt` in `preprocessing_v8`

### Select tissues and individuals for downstream analyses
```
Rscript preprocessing/filter_tissues_individuals.R
```
Must be run from the upper level directory of the repo (e.g., the location of this readme).

Generates `preprocessing_v8/gtex_2017-06-05_v8_design_passed.txt`, `preprocessing_v8/gtex_2017-06-05_v8_individuals_passed.txt`, `preprocessing_v8/gtex_2017-06-05_v8_tissues_passed.txt`, and `preprocessing_v8/gtex_2017-06-05_v8_normalized_expression.subset.txt.gz`. Also produces summary figures `figures/gtex_v8_design.pdf`. The subset file filtered for missingness is used in correlation-outlier calling. No missingness filter is applied for multi-tissue eOutlier calling.

