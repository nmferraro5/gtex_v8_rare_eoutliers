library(data.table)
library(dplyr)

gtex_dir = '/srv/scratch/restricted/GOATs/features_v8/byGene/10kb_genebody/'
data_dir = '/users/nferraro/data/goats_data/v8_data/'

cgenes = fread(paste0(data_dir, 'coloc/medz_sig_coloc_RV_genes.txt'),header=F)$V1
goutliers = fread(paste0(data_dir, 'gtexV8_global_outliers_medz3_iqr.txt'),header=F)$V1

rare_snps = fread(paste0(gtex_dir, 'all_rare_variants_SNPs.txt'),header=F) %>% 
  filter(!(V1 %in% goutliers)) %>% filter(V2 %in% cgenes) %>% mutate(Type='SNP')
rare_indels = fread(paste0(gtex_dir, 'all_rare_variants_indels.txt'),header=F) %>% 
  filter(!(V1 %in% goutliers)) %>% filter(V2 %in% cgenes) %>% mutate(Type='indel')
rare_SVs = fread(paste0(gtex_dir, 'all_rare_variants_HallLabSV.txt'),header=F) %>% 
  filter(!(V1 %in% goutliers)) %>% filter(V2 %in% cgenes) %>% mutate(Type='SV')

all_vars = rbind(rare_snps,rare_indels,rare_SVs)
colnames(all_vars) = c('Ind','Gene','Chr','Start','Stop','Type')

write.table(all_vars, file=paste0(data_dir, 'gtex_medz_coloc_genes_all_rare_vars.txt'),sep='\t',quote=F,row.names=F)

