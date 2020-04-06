Outlier calling
====================


For GTEx v8 samples
-------------------

### Set up directories for results
```
mkdir $RAREDIR/data_v8/outliers
```

### Run outlier calling
```
Rscript outlier_calling/call_outliers.R \
        --Z.SCORES=$RAREDIR/preprocessing_v8/gtex_2017-06-05_normalized_expression.txt.gz \
        --GLOBAL=${datadir}gtexV8_global_outliers_medz3_iqr.txt \
        --OUT.PREFIX=/data_v8/outliers/gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed \
        --N.PHEN=5 
```

Generates one file in specified directory with columns gene ind N Df MedZ Y, where N refers to the number of individuals tested for a given gene, Df refers to the number of measurements available for that gene in that individual and Y indicates outlier or control. The `N.PHEN` argument specifies the minimum number of measurements (i.e. tissues) available for a gene-individual required for inclusion. We use 5. `GLOBAL` specifies a file with a list of individual IDs to remove as global outliers. If NA, no individuals are removed.

### Identify global outlier individuals 
```
Rscript outlier_calling/identify_global_outliers.R \
        --OUTLIERS=$RAREDIR/data_v8/outliers/gtexV8.outlier.control.medz.txt \
        --METHOD=proportion
```
Removes global outlier individuals, either defined based on the proportion of genes called as outliers per individual relative to the population or the number of genes called as outliers, determined by setting 'proportion' or 'number' in `METHOD`. Writes out a new outlier file with `_globalOutliersRemoved` appended to the specified outlier txt file.

### From the multi-tissue outliers, determine which tissues have extreme effects and outlier sharing across tissues
```
Rscript outlier_calling/extract_extreme_tissues.R \
    --Z.SCORES=$RAREDIR/preprocessing_v8/gtex_2017-06-05_normalized_expression.txt.gz \
    --EXP.DATA=$RAREDIR//preprocessing_v8/gtex_2017-06-05_normalized_expression_v8ciseQTLs_removed.txt.gz \
```

Generates figures in `figures/GTEXv8_pair_jaccard.pdf`.
