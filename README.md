GTEx v8 eOutlier Analysis 
=====================================

Setup
-----
First set paths in `paths`.  
Then run to access the path variables:
```
source <path to repo>/paths
```

The setup script makes the basic directory structure:
```
bash 0_setup.sh
```
This structure is described below.
```
RAREDIR/
	data_v8/
	features_v8/
	figures/
	paper_figures/
	preprocessing_v8/
```

Other subdirectories are created in individual scripts later.

Data processing and Analysis
----------------------------
[Data correction, transformation, and reorganization](correction.md)  
[Feature construction](features.md)  
[Imputation and outlier calling](outlier_calling.md)  
[Rare variant enrichments](enrichment.md)  
