#!/bin/env/Rscript

library(data.table)
library(dplyr)

data_dir = '/users/nferraro/data/goats_data/v8_data/'
baseDir = Sys.getenv('RAREDIR')
print(baseDir)

sp_ws = fread(paste0(baseDir, '/data_v8/watershed/splicing.gam.watershed.rounded.variants.txt'))
exp_ws = fread(paste0(baseDir, '/data_v8/watershed/expression.gam.watershed.rounded.variants.txt'))
ase_ws = fread(paste0(baseDir, '/data_v8/watershed/ase.gam.watershed.rounded.variants.txt'))
ws_all = rbind(exp_ws %>% mutate(outlier_pvalue = total_expression_outlier_pvalue, gam_posterior = total_expression_gam_posterior, ws_posterior = total_expression_watershed_posterior) %>% 
                 select(-total_expression_outlier_pvalue, -total_expression_gam_posterior, -total_expression_watershed_posterior) %>% mutate(Type = 'TE'),
               sp_ws %>% mutate(outlier_pvalue = splicing_outlier_pvalue, gam_posterior = splicing_gam_posterior, ws_posterior = splicing_watershed_posterior) %>% 
                 select(-splicing_outlier_pvalue, -splicing_gam_posterior, -splicing_watershed_posterior) %>% mutate(Type = 'Splicing'),
               ase_ws %>% mutate(outlier_pvalue = ase_outlier_pvalue, gam_posterior = ase_gam_posterior, ws_posterior = ase_watershed_posterior) %>% 
                 select(-ase_outlier_pvalue, -ase_gam_posterior, -ase_watershed_posterior) %>% mutate(Type = 'ASE'))

print('Watershed')
ws_top = ws_all %>% group_by(VID,Type) %>% top_n(1,ws_posterior) %>% ungroup() %>% mutate(Model = 'Watershed')
print('GAM')
gam_top = ws_all %>% group_by(VID,Type) %>% top_n(1,gam_posterior) %>% ungroup() %>% mutate(Model = 'GAM')
out_top = rbind(ws_top, gam_top)
out_top$outlier_pvalue = round(out_top$outlier_pvalue, 4)
out_top$gam_posterior = round(out_top$gam_posterior, 4)
out_top$ws_posterior = round(out_top$ws_posterior, 4)
write.table(out_top, file=paste0(data_dir, 'coloc/gtexV8.top.watershed.gam.per.type.per.variant.txt'), sep='\t', quote=F, row.names=F)

#eqtl_coloc_results = fread(paste0(data_dir, 'coloc/coloc_sig_results_0.5.txt'))
#eqtl_coloc_results$Trait = sapply(eqtl_coloc_results$filename, function(x) strsplit(x, '__PM__')[[1]][1])
#
#eqtl_results = fread(paste0(data_dir, 'coloc/enloc_sig_results_0.5.txt'))
#eqtl_results$Trait = sapply(eqtl_results$filename, function(x) strsplit(x, '_w_')[[1]][2])
#eqtl_results$Trait = sapply(eqtl_results$Trait, function(x) strsplit(x, '_enloc_output.txt')[[1]][1])
#sqtl_results = fread(paste0(data_dir, 'splicing/sqtl_enloc_rcp0.5_collapsed_bygene.txt'))
#colnames(sqtl_results) = c('Chr', 'Start', 'End', 'rcp', 'Trait', 'Tissue', 'GeneChr', 'GeneStart', 'GeneEnd', 'Gene')
#sqtl_results = sqtl_results %>% select(rcp, Trait, Gene)
#
#eqtl_results = eqtl_results %>% select(rcp, Trait, gene_id)
#colnames(eqtl_results) = c('ColocProb', 'Enloc_name', 'Gene')
#eqtl_coloc_results = eqtl_coloc_results %>% select(p4,Trait,gene)
#colnames(eqtl_coloc_results) = c('ColocProb', 'Coloc_name', 'Gene')
#colnames(sqtl_results) = c('ColocProb', 'Sqtl_name', 'Gene')
#
#coloc_genes = unique(c(eqtl_coloc_results$Gene, eqtl_results$Gene, sqtl_results$Gene))

