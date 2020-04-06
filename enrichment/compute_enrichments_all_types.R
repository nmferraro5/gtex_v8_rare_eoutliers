#!/usr/bin/env Rscript

# Single tissue outliers

## Get path to the data directories
dir = Sys.getenv('RAREDIR')

## Load required packages
require(data.table)
require(dplyr)
library(ggplot2)
library(epitools)
library(RColorBrewer)
library(argparse)

source('enrichment/enrichment_functions_MEDZ.R')

parser = ArgumentParser()
parser$add_argument('--pvalue', help = 'Outlier threshold')
args = parser$parse_args()

pt = as.numeric(args$pvalue)

data_dir = '/users/nferraro/data/goats_data/v8_data/'
window = '10kb_genebody'
featdir = paste0(dir, '/features_v8/byGene/', window, '/')

medz_data = fread(paste0(data_dir, 'gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.medz.filtered.txt'))
medz_data$Pvalue = 2*pnorm(-abs(medz_data$MedZ))

#medz_outliers = filter(medz_data, abs(MedZ) > zt) %>% mutate(Y=1)
#medz_controls = filter(medz_data, abs(MedZ) < zt, Gene %in% medz_outliers$Gene) %>% mutate(Y=0)
medz_outliers = filter(medz_data, Pvalue < pt) %>% mutate(Y=1)
print(nrow(medz_outliers))
medz_controls = filter(medz_data, Pvalue > pt, Gene %in% medz_outliers$Gene) %>% mutate(Y=0)
print(nrow(medz_controls))
medz_data = rbind(medz_outliers, medz_controls) %>% select(-Pvalue)
print(nrow(medz_data))
colnames(medz_data)[3] = 'OutlierValue'
medz_data$Gene = sapply(medz_data$Gene, function(x) strsplit(x,'[.]')[[1]][1])
medz_data = medz_data %>% mutate(UID = paste(Ind,Gene,sep='_')) %>% mutate(Method = 'TotalExpression')

splice_data = fread(paste0('zcat ', data_dir, 'splicing/cross_tissue_covariate_method_none_no_global_outliers_ea_only_emperical_pvalue_gene_level.txt.gz'))
splice_melted = melt(splice_data)
rm(splice_data)
splice_outliers = filter(splice_melted, value < pt) %>% mutate(Y=1)
splice_controls = filter(splice_melted, value > pt, CLUSTER_ID %in% splice_outliers$CLUSTER_ID) %>% mutate(Y=0)
splice_data = rbind(splice_outliers,splice_controls)
rm(splice_melted)
colnames(splice_data)[1:3] = c('Gene', 'Ind', 'OutlierValue')
splice_data$Gene = sapply(splice_data$Gene, function(x) strsplit(x, '[.]')[[1]][1])
splice_data = filter(splice_data, !is.na(OutlierValue))
splice_data = splice_data %>% mutate(UID = paste(Ind,Gene,sep='_')) %>% mutate(Method = 'Splicing')

ase_data = fread(paste0('zcat ', data_dir, 'ASE/median.ad.scores.uncorrected.no.global.outliers.v8.vg.tsv.gz'))
ase_melted = melt(ase_data)
rm(ase_data)
colnames(ase_melted)[1] = 'GeneID'
ase_outliers = filter(ase_melted, value < pt) %>% mutate(Y=1)
ase_controls = filter(ase_melted, value > pt, GeneID %in% ase_outliers$GeneID) %>% mutate(Y=0)
ase_data = rbind(ase_outliers,ase_controls)
#ase_data = ase_melted
rm(ase_melted)
colnames(ase_data)[1:3] = c('Gene','Ind', 'OutlierValue')
ase_data = filter(ase_data, !is.na(OutlierValue))
ase_data = ase_data %>% mutate(UID = paste(Gene,Ind,sep='_')) %>% mutate(Method = 'ASE')

outlier_list = list(medz_data, splice_data, ase_data)

all_logits = do.call(rbind, lapply(outlier_list, get.all.enrich, featdir, counts = TRUE))

write.table(all_logits, file=paste0(data_dir, 'gtexV8.all.data.types.svs.thresholds.', as.character(pt), '.10kb.gnomad.txt'),sep='\t',quote=F,row.names=F)


