#!/usr/bin/env Rscript

rm(list = ls())

library(data.table)
library(dplyr)
library(stringr)
library(argparse)

## Script for combining outlier calls and features then computing enrichments.

baseDir = Sys.getenv('RAREDIR')

source('enrichment/enrichment_functions.R')

parser = ArgumentParser()
parser$add_argument('--dir.suffix', help = 'Suffix of directory name')
parser$add_argument('--outliers.prefix', help = 'Prefix of outlier files.')
parser$add_argument('--window', help = 'Window name.')
parser$add_argument('--fdr.thresh', default = 0.01, help = 'FDR threshold.')
parser$add_argument('--z.thresh', default = 2, help = 'Z-score threshold.')
parser$add_argument('--output.suffix', default = '', help = 'suffix for output file')
parser$add_argument('--eqtl.correct', default = 'False', help = 'Use ciseQTL corrected file')
args = parser$parse_args()

iodir = paste0(baseDir, '/data_', args$dir.suffix, '/outliers/')
featdir = paste0(baseDir, '/features_', args$dir.suffix, '/byGene/', args$window, '/')
## Get files with outliers
outlier_files = list.files(iodir, paste0(args$outliers.prefix, '.+txt'), full.names = T)
print(outlier_files)
#kind = which(grepl('zero',outlier_files) | grepl('none',outlier_files))
#print(outlier_files[kind])
#outlier_files = outlier_files[kind]
#outlier_files = outlier_files[c(2,4,8)]
if (args$eqtl.correct == 'False') {
  kinds = which(grepl('ciseQTLs', outlier_files))
  outlier_files = outlier_files[-kinds]
}

#outlier_files = c('/users/nferraro/data/goats_data/outliers_norm_data/post_eqtl_correction/outliers.rp.zscore.knn.txt',
#			'/users/nferraro/data/goats_data/outliers_norm_data/post_eqtl_correction/outliers.rp.zscore.medz.txt')

# Restricting just to KNN, MEDZ, and STFZ
#rinds = c(grep('soft[.]txt', outlier_files), grep('em[.]txt', outlier_files), grep('mean[.]txt', outlier_files), grep('pmd[.]txt', outlier_files))
#outlier_files = outlier_files[-rinds]

#medz_indx = grep('medz', outlier_files)
args$z.thresh = as.numeric(args$z.thresh)

outlier_file = paste0(baseDir, '/data_v8/outliers/gtexV8.outlier.controls.medz.txt')
## First get medz outliers
#medz_outliers = filter.outliers(outlier_files[medz_indx], args$z.thresh)
medz_outliers = filter.outliers(outlier_file, args$z.thresh)
noutliers = sum(medz_outliers$Y == 'outlier')
print(noutliers)
#noutliers = 3281
outlier_list_top = list(medz_outliers)

## then outliers from other methods (either by FDR or matching the number of medz outliers) and add medz ones
#outlier_list_fdr = lapply(outlier_files[-medz_indx], filter.outliers, args$fdr.thresh)
#outlier_list_fdr[[length(outlier_list_fdr) + 1]] = medz_outliers

#outlier_list_top = lapply(outlier_files[-medz_indx], filter.outliers, noutliers, topn = TRUE)
#outlier_list_top = lapply(outlier_files, filter.outliers, noutliers, topn = TRUE)
#outlier_list_top[[length(outlier_list_top) + 1]] = medz_outliers
# Get log odds ratio data for all methods (both features and counts)
#all_logits_fdr = rbind(do.call(rbind, lapply(outlier_list_fdr, get.all.enrich, featdir)),
#                       do.call(rbind, lapply(outlier_list_fdr, get.all.enrich, featdir, counts = TRUE)))

#all_logits_top = rbind(do.call(rbind, lapply(outlier_list_top, get.all.enrich, featdir)),
#                        do.call(rbind, lapply(outlier_list_top, get.all.enrich, featdir, counts = TRUE)))

all_logits_top = do.call(rbind, lapply(outlier_list_top, get.all.enrich, featdir, counts = TRUE))

#all_logits_fdr = do.call(rbind, lapply(outlier_list_fdr, get.all.enrich, featdir, counts = TRUE))


### BREAK
#all_logits_top = do.call(rbind, lapply(outlier_list_top, get.all.enrich, featdir, counts = TRUE))

## Set factor levels and labels for plotting
#maf.levels = unique(all_logits_fdr$Maf)
#maf.labels = sub('MAF', '', maf.levels)
#maf.order = order(as.numeric(str_split_fixed(maf.labels, '-', 2)[,1])) # get order by first number
#maf.levels = maf.levels[maf.order]
#maf.labels = maf.labels[maf.order]
#all_logits_fdr$Maf = factor(all_logits_fdr$Maf, levels = maf.levels, labels = maf.labels)
#all_logits_top$Maf = factor(all_logits_top$Maf, levels = maf.levels, labels = maf.labels)

## Save relevant objects
#save(outlier_list_fdr, file = 'gtexV8.fdr.outliers.cor.RData')
save(outlier_list_top, all_logits_top,
     file = paste0(iodir, 'enrichments_relativeRisk_', args$window, '_Z', args$z.thresh, '_FDR', args$fdr.thresh, args$output.suffix, '.RData'))
