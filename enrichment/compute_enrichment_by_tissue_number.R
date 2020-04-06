#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(epitools)

goats_dir = Sys.getenv('RAREDIR')
data_dir = '/users/nferraro/data/goats_data/v8_data/'

exp_outliers = fread(paste0(goats_dir, '/data_v8/outliers/gtexV8.expression.outliers.across.tissue.numbers.N5.txt')) %>% select(Gene,Tissue,variable,value,NT,NO)
ase_outliers = fread(paste0(goats_dir, '/data_v8/ase_v8/gtexV8.ase.outliers.across.tissue.numbers.N5.txt'))
gene_names = fread(paste0(data_dir, 'gencode.v26.GRCh38.genes.bed')) %>% select(V4)
gene_names$GeneID = sapply(gene_names$V4, function(x) strsplit(x, '[.]')[[1]][1])
colnames(gene_names)[1] = 'Gene'
colnames(ase_outliers)[1] = 'GeneID'
ase_outliers = inner_join(ase_outliers, gene_names, by='GeneID') %>% select(Gene,Tissue,variable,value,NT,NO)
sp_outliers = fread(paste0(goats_dir, '/data_v8/splicing/gtexV8.splicing.outliers.across.tissue.numbers.N5.txt'))

rare_snps = fread(paste0(data_dir, 'all_genes_with_rare_SNPs.txt'), header=F)
rare_indels = fread(paste0(data_dir, 'all_genes_with_rare_indels.txt'), header=F)
rare_si = rbind(rare_snps, rare_indels) %>% distinct() 
si_inds = unique(rare_si$V1)
rm(rare_snps, rare_indels)
rare_svs = fread(paste0(data_dir, 'all_genes_with_rare_HallLabSV.txt'), header=F)
rare_si_uids = paste(rare_si$V1, rare_si$V2, sep=':')
rare_sv_uids = paste(rare_svs$V1, rare_svs$V2, sep=':')
sv_inds = unique(rare_svs$V1)
rm(rare_si, rare_svs)

get_enrichments <- function(tt, data, dcat) {
  data = filter(data, variable %in% si_inds, NT >= 5, !is.na(value)) %>% distinct()
  if (tt == 1 | tt == 2) {
    outliers = filter(data, NO >= tt) %>% mutate(UID = paste(variable,Gene,sep=':')) %>% 
      mutate(HasSI = ifelse(UID %in% rare_si_uids, 1, 0), HasSV = ifelse(UID %in% rare_sv_uids, 1, 0))
    controls = filter(data, NO < tt, Gene %in% outliers$Gene) %>% mutate(UID = paste(variable,Gene,sep=':')) %>% 
      mutate(HasSI = ifelse(UID %in% rare_si_uids, 1, 0), HasSV = ifelse(UID %in% rare_sv_uids, 1, 0))

  } else {
    outliers = filter(data, NO >= (tt/100)*NT) %>% mutate(UID = paste(variable,Gene,sep=':')) %>% 
      mutate(HasSI = ifelse(UID %in% rare_si_uids, 1, 0), HasSV = ifelse(UID %in% rare_sv_uids, 1, 0))
    controls = filter(data, NO < (tt/100)*NT, Gene %in% outliers$Gene) %>% mutate(UID = paste(variable,Gene,sep=':')) %>% 
      mutate(HasSI = ifelse(UID %in% rare_si_uids, 1, 0), HasSV = ifelse(UID %in% rare_sv_uids, 1, 0))
  }
  si_table = rbind(c(nrow(filter(controls, HasSI == 0)), nrow(filter(controls, HasSI == 1))),
                     c(nrow(filter(outliers, HasSI == 0)), nrow(filter(outliers, HasSI == 1))))
  sv_table = rbind(c(nrow(filter(controls, variable %in% sv_inds, HasSV == 0)), nrow(filter(controls, variable %in% sv_inds, HasSV == 1))),
                     c(nrow(filter(outliers, variable %in% sv_inds, HasSV == 0)), nrow(filter(outliers, variable %in% sv_inds, HasSV == 1))))
  si_rr = epitab(si_table, method='riskratio')$tab
  sv_rr = epitab(sv_table, method='riskratio')$tab
  estims = data.frame(Riskratio = si_rr[2,5], Lower = si_rr[2,6], Upper = si_rr[2,7], Pval = si_rr[2,8], Feature = 'SNVs+indels', Threshold = tt, Category = dcat)
  estims = rbind(estims, data.frame(Riskratio = sv_rr[2,5], Lower = sv_rr[2,6], Upper = sv_rr[2,7], Pval = sv_rr[2,8], Feature = 'SVs', Threshold = tt, Category = dcat))
  return(estims)
}

thresholds = c(1, 2, 25, 50, 75, 90, 100)
ase_estims = do.call(rbind, lapply(thresholds, function(x) get_enrichments(x, ase_outliers, 'aseOutliers')))
exp_estims = do.call(rbind, lapply(thresholds, function(x) get_enrichments(x, exp_outliers, 'eOutliers')))
sp_estims = do.call(rbind, lapply(thresholds, function(x) get_enrichments(x, sp_outliers, 'sOutliers')))

all_estims = rbind(ase_estims, exp_estims, sp_estims)
write.table(all_estims, file=paste0(data_dir, 'enrichments/gtexV8.rare.enrichments.by.number.tissues.all.data.types.txt'), sep='\t', quote=F, row.names=F)



