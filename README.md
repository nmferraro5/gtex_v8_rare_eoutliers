GTEx v8 Outlier Analysis 
=====================================

Setup
-----
First set paths in `paths`.  
Then run to access the path variables:
```
source <path to repo>/paths
```

Then you can run the setup script with makes the basic directory structure for you:
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
