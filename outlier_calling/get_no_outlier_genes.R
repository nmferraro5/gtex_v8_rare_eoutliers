library(data.table)
library(dplyr)

goats_dir = '/srv/scratch/restricted/GOATs/'
data_dir = '/users/nferraro/data/goats_data/v8_data/'

te_files = dir(paste0(goats_dir, 'preprocessing_v8/PEER_v8/'), '*.peer.v8ciseQTLs.ztrans.txt')
ase_files = dir(paste0(goats_dir, 'data_v8/ase_v8/'), 'combined.uncorrected.no.global.outlier.genes.no.global.outlier.indvs.ad.scores.in.*.tsv.gz')
sp_files = dir(paste0(goats_dir, 'data_v8/splicing/'), '*_covariate_method_none_no_global_outliers_ea_only_merged_emperical_pvalue_gene_level.txt')

all_genes = all_genes = fread(paste0(data_dir, 'gencode.v26.GRCh38.genes_genetypes_autosomal_PCandlinc_only.txt'),header=F)

te_no_outlier_genes = all_genes$V1

for (tf in te_files) {
  te_data = melt(fread(paste0(goats_dir, 'preprocessing_v8/PEER_v8/', tf)))
  outlier_genes = unique(filter(te_data, abs(value) > 3)$Id)
  kinds = which(!(te_no_outlier_genes %in% outlier_genes))
  te_no_outlier_genes = te_no_outlier_genes[kinds]
}

ase_no_outlier_genes = sapply(all_genes$V1, function(x) strsplit(x, '[.]')[[1]][1])

for (af in ase_files) {
  ase_data = melt(fread(paste0('zcat ', goats_dir, 'data_v8/ase_v8/', af)))
  outlier_genes = unique(filter(ase_data, value < 0.0027)$GeneID)
  kinds = which(!(ase_no_outlier_genes %in% outlier_genes))
  ase_no_outlier_genes = ase_no_outlier_genes[kinds]
}

sp_no_outlier_genes = all_genes$V1

for (sf in sp_files) {
  sp_data = melt(fread(paste0(goats_dir, 'data_v8/splicing/', sf)))
  outlier_genes = unique(filter(sp_data, value < 0.0027)$GENE_ID)
  kinds = which(!(sp_no_outlier_genes %in% outlier_genes))
  sp_no_outlier_genes = sp_no_outlier_genes[kinds]
}

save(te_no_outlier_genes, ase_no_outlier_genes, sp_no_outlier_genes, file=paste0(data_dir, 'example_genes/no_outlier_genes.RData'))

