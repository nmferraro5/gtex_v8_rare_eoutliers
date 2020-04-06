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
parser$add_argument('--dir.suffix', help = 'Suffix of directory name (e.g., v7 or multi_omics).')
parser$add_argument('--outliers.file', help = 'Outlier file.')
parser$add_argument('--window', help = 'Window name.')
parser$add_argument('--z.thresh', default = 3, help = 'Z-score threshold.')
parser$add_argument('--output.suffix', default = '', help = 'suffix for output file')
args = parser$parse_args()

iodir = paste0(baseDir, '/data_', args$dir.suffix, '/outliers/')
featdir = paste0(baseDir, '/features_', args$dir.suffix, '/byGene/', args$window, '/')

## Get files with outliers
args$z.thresh = as.numeric(args$z.thresh)

## First get medz outliers
medz_outliers = fread(paste0(baseDir, '/data_v8/outliers/', args$outliers.file)) %>% mutate(Method = 'MEDZ')

all_logits_top = get.all.enrich(medz_outliers, featdir, counts = TRUE)

## Set factor levels and labels for plotting
maf.levels = unique(all_logits_top$Maf)
maf.labels = sub('MAF', '', maf.levels)
maf.order = order(as.numeric(str_split_fixed(maf.labels, '-', 2)[,1])) # get order by first number
maf.levels = maf.levels[maf.order]
maf.labels = maf.labels[maf.order]
all_logits_top$Maf = factor(all_logits_top$Maf, levels = maf.levels, labels = maf.labels)

save(all_logits_top, file = paste0(iodir, 'enrichments_', args$window, '_Z', args$z.thresh, '_', args$output.suffix, '.RData'))

