#!/usr/bin/env

## Merge outliers with variant annotations 

library(data.table)
library(dplyr)
library(ggplot2)
require(RColorBrewer)
library(argparse)

RAREDIR = Sys.getenv('RAREDIR')

parser = ArgumentParser()
parser$add_argument('--method', help = 'Method')
args = parser$parse_args()

method = args$method

get_bin <- function(pval) {
  if (pval < 1e-13) {
    bin = '0~1e-13'
  } else if (pval < 1e-11) {
    bin = '1e-13~1e-11'
  } else if (pval < 1e-09) {
    bin = '1e-11~1e-09'
  } else if (pval < 1e-07) {
    bin = '1e-09~1e-07'
  } else if (pval < 1e-05) {
    bin = '1e-07~1e-05'
  } else if (pval < 1e-04) {
    bin = '1e-05~1e-04'
  } else if (pval < 1e-03) {
    bin = '1e-04~1e-03'
  } else if (pval < 1e-02) {
    bin = '1e-03~1e-02'
  } else if (pval < 5e-02) {
    bin = '1e-02~5e-02'
  } else {
    bin = 'nonOutlier'
  }
  return(bin)
}

if (method == 'medz') {
  medz_outliers = fread(paste0(RAREDIR, '/data_v8/outliers/gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.medz.txt'))
  colnames(medz_outliers)[1:2] = c('indiv_id','gene_id')
  medz_outliers$indiv_id = sapply(medz_outliers$indiv_id, function(x) strsplit(x,'GTEX-')[[1]][2])
} else if (method == 'splicing') {
  splicing_outliers = fread(paste0(RAREDIR, '/data_v8/splicing/cross_tissue_covariate_method_none_no_global_outliers_ea_only_emperical_pvalue_gene_level.txt'))
  splicing_melted = melt(splicing_outliers)
  splicing_melted = filter(splicing_melted,!is.nan(value))
  splicing_bins = sapply(splicing_melted$value, function(x) get_bin(x))
  splicing_melted$pval_bin = splicing_bins
  colnames(splicing_melted)[1:2] = c('gene_id','indiv_id')
  splicing_melted$indiv_id = as.character(splicing_melted$indiv_id)
  splicing_melted$indiv_id = sapply(splicing_melted$indiv_id, function(x) strsplit(x,'GTEX-')[[1]][2])
  outlier_genes = filter(splicing_melted,pval_bin!='nonOutlier')$gene_id
  splicing_melted = filter(splicing_melted,gene_id %in% outlier_genes)
} else if (method == 'ase') {
  #ase_data = fread(paste0(RAREDIR, '/data_v8/ase/combined.ad.scores.in.MEDIAN_4_10_update.tsv'))
  ase_data = fread('/users/nferraro/data/goats_data/v8_data/ASE/median.ad.scores.uncorrected.no.global.outliers.v8.vg.tsv', stringsAsFactors=F)
  ase_outliers = melt(ase_data)
  ase_outliers$pval_bin = sapply(ase_outliers$value, function(x) ifelse(is.na(x), NA, get_bin(x)))
  colnames(ase_outliers)[1:2] = c('gene_id','indiv_id')
  ase_outliers$indiv_id = as.character(ase_outliers$indiv_id)
  ase_outliers$indiv_id = sapply(ase_outliers$indiv_id, function(x) strsplit(x,'GTEX-')[[1]][2])
  gene_map = fread('/users/nferraro/data/goats_data/v8_data/gencode.v26.GRCh38.genes.bed',header=F)
  gene_map$gene_id = sapply(gene_map$V4, function(x) strsplit(x, '[.]')[[1]][1])
  gene_map = gene_map %>% select(V4,gene_id)
  colnames(gene_map) = c('Gene', 'gene_id')
  ase_outliers = inner_join(ase_outliers, gene_map, by='gene_id')
  ase_outliers = ase_outliers %>% select(Gene,indiv_id,value,pval_bin)
  colnames(ase_outliers)[1] = 'gene_id'
  
} else {
  cor_data = fread(paste0(RAREDIR, '/data_v8/outliers/gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.subset.newParams.knn.txt'))
  cor_outliers = filter(cor_data,Y=='outlier',FDR < 0.01)
  cor_controls = filter(cor_data, Gene %in% cor_outliers$Gene, Y=='control')
  cor_bins = sapply(cor_outliers$FDR, function(x) get_bin(x))
  cor_outliers$pval_bin = cor_bins
  cor_controls$pval_bin = 'nonOutlier'
  cor_outliers = rbind(cor_outliers %>% select(Ind,Gene,FDR,pval_bin,Y),cor_controls %>% select(Ind,Gene,FDR,pval_bin,Y))
  colnames(cor_outliers)[1:2] = c('indiv_id','gene_id')
  cor_outliers$indiv_id = sapply(cor_outliers$indiv_id, function(x) strsplit(x,'GTEX-')[[1]][2])
}

variant_file = read.table('/users/xli6/projects/gtex/annotation/combined/gtex_v8_rare_GxI_collapsed_feature.tsv',header=T,fill=NA)
variant_file$gene_id = as.character(variant_file$gene_id)
if (method == 'medz') {
  medz_variants = merge(medz_outliers,variant_file,by=c('indiv_id','gene_id'),all.x=T)
  write.table(medz_variants,file=paste0(RAREDIR,'/data_v8/outliers/gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.medz.lncRNA.variantAnnotations.withPromoter.txt'),sep='\t',quote=F,row.names=F)
} else if (method == 'splicing') {
  splice_variants = merge(splicing_melted,variant_file,by=c('indiv_id','gene_id'),all.x=T)
  write.table(splice_variants,file=paste0(RAREDIR,'/data_v8/splicing/gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.splice.variantAnnotations.withPromoter.txt'),sep='\t',quote=F,row.names=F)
} else if (method == 'ase') {
  print(ase_outliers[1,])
  print(variant_file[1,])
  ase_variants = merge(ase_outliers,variant_file,by=c('indiv_id','gene_id'),all.x=T)
  write.table(ase_variants,file=paste0(RAREDIR,'/data_v8/ase/gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.ase.variantAnnotations.v8.vg.txt'),sep='\t',quote=F,row.names=F)
} else {
  cor_variants = merge(cor_outliers,variant_file,by=c('indiv_id','gene_id'),all.x=T)
  write.table(cor_variants,file=paste0(RAREDIR,'/data_v8/outliers/gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.knn.variantAnnotations.txt'),sep='\t',quote=F,row.names=F)
}

