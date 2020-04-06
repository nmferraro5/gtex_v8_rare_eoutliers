#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(ggplot2)
require(RColorBrewer)

### sqtl analysis
data_dir = '/users/nferraro/data/goats_data/v8_data/splicing/'
sqtl_results = fread(paste0(data_dir, 'predixcan_splicing_significant.csv'))
sqtl_results$Start = sapply(sqtl_results$gene, function(x) strsplit(x,'_')[[1]][3])
sqtl_results$End = sapply(sqtl_results$gene, function(x) strsplit(x,'_')[[1]][4])
sqtl_results$Start = as.numeric(sqtl_results$Start)
sqtl_results$End = as.numeric(sqtl_results$End)
gene_positions = read.table('/users/nferraro/data/goats_data/v8_data/gencode.v26.GRCh38.genes.gtf',fill=NA)
gene_positions = gene_positions %>% filter(V3 == 'gene') %>% select(V1,V4,V5,V10)

splicing_variants = fread(paste0(data_dir, 'gtexV8.outliers.splicing.withRareVars.txt'))
splicing_outliers = filter(splicing_variants,value < 0.0027)
sqtl_results$chr = paste0('chr',sqtl_results$chr)

splicing_genes = filter(gene_positions,V10 %in% splicing_outliers$gene_id)
window = 10000
test_in_window <- function(gstart, gend, istart, iend,window) {
  ## test within gene
  if (istart > gstart - window && iend < gend + window) {
    return(1)
  } else {
    return(0)
  }
}

get_sqtl_gene_overlap <- function(i) {
  gstart = as.numeric(as.character(splicing_genes$V4[i]))
  gend = as.numeric(as.character(splicing_genes$V5[i]))
  gchr = as.character(splicing_genes$V1[i])
  sqtl_chr = filter(sqtl_results, chr == gchr)
  kinds = sapply(1:nrow(sqtl_chr), function(x) ifelse(sqtl_chr$pos1[x] > gstart && sqtl_chr$pos2[x] < gend, 1, 0))
  gene_sqtls = sqtl_chr[which(kinds==1),]
  gene_sqtls$gene_id = rep(as.character(splicing_genes$V10[i]),nrow(gene_sqtls))
  gene_sqtls$GStart = rep(gstart,nrow(gene_sqtls))
  gene_sqtls$GEnd = rep(gend,nrow(gene_sqtls))
  return(gene_sqtls)
}

sqtl_gene_overlap = do.call(rbind, lapply(1:nrow(splicing_genes),get_sqtl_gene_overlap))
splicing_outliers = merge(splicing_outliers,sqtl_gene_overlap,by=c('gene_id'),all.x=T)



