Feature construction
====================

All paths for created files are relative to `$RAREDIR` unless otherwise noted.

Basic processing of annotation files
------------------------------------
```
bash preprocessing_v8/process_reference_files.sh $GTEX_GENES
```
Creates padded bed files `preprocessing_v8/gencode.v26.GRCh38.genes_padded10kb.bed` and `preprocessing_v8/gencode.v26.GRCh38.genes_padded10kb_PCandlinc_only.bed` and a filtered autosomal file `preprocessing_v8/gencode.v26.GRCh38.genes_genetypes_autosomal.txt` and `preprocessing_v8/gencode.v26.GRCh38.genes_genetypes_autosomal_PCandlinc_only.txt`.

For GTEx v8 samples
-------------------
###Make bed files for each individual from VCF including variants with MAF<=0.25 (separate ones for SNPs and indels). Also takes care of adding CADD scores to SNP ones.
```
bash feature_construction/vcf2bedfiles.sh &> feature_construction/vcf2bedfiles.log
```
relies on `vcf2bedfiles_processVCF.sh`, `vcf2bedfiles_processVCFtoolsOutput.sh`,
`vcf2bedfiles_compile_CADD_scores.sh`, and `vcf2bedfiles_compile_CADD_scores.py`

This creates the files under `features_v8/variantBeds` as well as
the lists of european individuals and the subset of them with WGS data:
`preprocessing/gtex_2016-01-15_v7_euro_ids.txt` and `preprocessing/gtex_2016-01-15_v7_euro_VCFids.txt`.

###Add CADD and other features to SNP and indel bed files.
Do this both for general features and individual tissue enhancer features.
```
bash feature_construction/run_add_features_variant_beds.sh &> feature_construction/run_add_features_variant_beds.log
bash feature_construction/run_add_enhancers_variant_beds.sh &> feature_construction/run_add_enhancers_variant_beds.log
```
Rely on `add_features_variant_beds.sh` and `add_enhancers_variant_beds.sh`.

Create files under `features_v8/bySite`.

###Collapse features by gene, grouping by allele frequency (both annotation features and counts).
```
bash feature_construction/run_build_feature_count_summaries_all_genes.sh &> feature_construction/run_build_feature_count_summaries_all_genes.log
bash feature_construction/run_build_enhancer_summaries_all_genes.sh &> feature_construction/run_build_enhancer_summaries_all_genes.log
```
Rely on `build_feature_summaries.py`, `build_feature_summaries_all_genes.sh`, `build_count_summaries_all_genes.sh`, `build_enhancer_summaries.py` and `build_enhancer_summaries_all_genes.sh`.
Can modify the scripts to generate features with different window sizes/variant filters.

Creates files under `features_v8/byGene`.

###Create bed files with features in exclusive distance bins
```
bash run_build_exclusive_count_summaries_all_genes.sh &> exclusive_summaries.log
```
Generates files under `features_v8/byGene` as above, using different exclusive regional windows upstream.


