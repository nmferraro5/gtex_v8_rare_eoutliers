#!/usr/bin/env Rscript

rm(list = ls())

library(data.table)
library(dplyr)
library(stringr)
library(argparse)

## Script for combining outlier calls and features then computing enrichments.

baseDir = Sys.getenv('RAREDIR')

source('enrichment/enrichment_functions_MEDZ_dist.R')

parser = ArgumentParser()
parser$add_argument('--dir.suffix', help = 'Suffix of directory name (e.g., v7 or multi_omics).')
parser$add_argument('--outliers.prefix', help = 'Prefix of outlier files.')
parser$add_argument('--window', help = 'Window name.')
parser$add_argument('--fdr.thresh', default = 0.01, help = 'FDR threshold.')
parser$add_argument('--z.thresh', default = 2, help = 'Z-score threshold.')
parser$add_argument('--output.suffix', default = '', help = 'suffix for output file')
args = parser$parse_args()

iodir = paste0(baseDir, '/data_', args$dir.suffix, '/outliers/')
featdir = paste0(baseDir, '/features_', args$dir.suffix, '/byGene/', args$window, '/')

## Get files with outliers
outlier_files = list.files(iodir, paste0(args$outliers.prefix, '.+txt'), full.names = T)
rinds = grep('subset',outlier_files)
outlier_files = outlier_files[-rinds]

# Restricting just to KNN, MEDZ, and STFZ
#rinds = c(grep('soft[.]txt', outlier_files), grep('em[.]txt', outlier_files), grep('mean[.]txt', outlier_files), grep('pmd[.]txt', outlier_files))
#outlier_files = outlier_files[-rinds]

medz_indx = grep('medz', outlier_files)
#medz_indx = grep('gtex.v6p.outliers.controls.medz', outlier_files)
zt = as.numeric(args$z.thresh)

data_dir = '/users/nferraro/data/goats_data/v8_data/'
## First get medz outliers
medz_data = fread(paste0(data_dir, 'gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.medz.filtered.txt'))

medz_outliers = filter(medz_data, abs(MedZ) > zt) %>% mutate(Y='outlier')
medz_controls = filter(medz_data, abs(MedZ) < zt, Gene %in% medz_outliers$Gene) %>% mutate(Y='control')

medz_data = rbind(medz_outliers, medz_controls)

medz_data = medz_data %>% mutate(Method = 'TotalExpression')

## ase outliers
pt = 2*pnorm(-abs(3))
ase_data = fread(paste0('zcat ', data_dir, 'ASE/median.ad.scores.uncorrected.no.global.outliers.v8.vg.tsv.gz'))
ase_melted = melt(ase_data)
colnames(ase_melted)[1] = 'GeneID'
rm(ase_data)
ase_outliers = filter(ase_melted, value < pt) %>% mutate(Y='outlier')
ase_controls = filter(ase_melted, value > pt, GeneID %in% ase_outliers$GeneID) %>% mutate(Y='control')
ase_data = rbind(ase_outliers,ase_controls)
rm(ase_melted)
colnames(ase_data)[1:2] = c('Gene','Ind')
ase_data = ase_data %>% mutate(Method = 'ASE')

### splicing outliers
splice_data = fread(paste0('zcat ', data_dir, 'splicing/cross_tissue_covariate_method_none_no_global_outliers_ea_only_emperical_pvalue_gene_level.txt.gz'))
splice_melted = melt(splice_data)
rm(splice_data)
splice_outliers = filter(splice_melted, value < pt) %>% mutate(Y='outlier')
splice_controls = filter(splice_melted, value > pt, CLUSTER_ID %in% splice_outliers$CLUSTER_ID) %>% mutate(Y='control')
splice_data = rbind(splice_outliers,splice_controls)
rm(splice_melted)
colnames(splice_data)[1:2] = c('Gene', 'Ind')
splice_data = splice_data %>% mutate(Method = 'Splicing')

#all_logits_top = get.all.enrich(medz_outliers, featdir, counts = TRUE)
outlier_list = list(medz_data, splice_data, ase_data)

all_logits_top = do.call(rbind, lapply(outlier_list, get.all.enrich, featdir, counts = TRUE))

## Set factor levels and labels for plotting
maf.levels = unique(all_logits_top$Maf)
maf.labels = sub('MAF', '', maf.levels)
maf.order = order(as.numeric(str_split_fixed(maf.labels, '-', 2)[,1])) # get order by first number
maf.levels = maf.levels[maf.order]
maf.labels = maf.labels[maf.order]
#all_logits_fdr$Maf = factor(all_logits_fdr$Maf, levels = maf.levels, labels = maf.labels)
all_logits_top$Maf = factor(all_logits_top$Maf, levels = maf.levels, labels = maf.labels)

## Save relevant objects
#save(outlier_list_fdr, outlier_list_top, all_logits_fdr, all_logits_top,
#     file = paste0(iodir, 'enrichments_', args$window, '_Z', args$z.thresh, '_FDR', args$fdr.thresh, args$output.suffix, '.RData'))
save(all_logits_top, file = paste0(iodir, 'enrichmentsDIST_downstream_relativeRisk_v8_vg_', args$window, '_Z', zt, '_FDR', args$fdr.thresh, args$output.suffix, '.RData'))
