#!/usr/bin/env Rscript

rm(list = ls())

library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)
require(foreach)
require(doMC)
require(RColorBrewer)
library(epitools)
#library(argparse)

registerDoMC(cores = 1)

## Script for combining outlier calls and features then computing enrichments.

baseDir = Sys.getenv('RAREDIR')

source('enrichment_functions.R')

###### FUNCTIONS
## Takes as input a data frame with all the tissue specific outliers for that method.
## Extracts outliers (and matched controls) for the given tissue
## Returns the relevant data frame subset
get.tissue.outliers.controls <- function(tissue, outliers.controls, zt=NA) {
  outliers.subset = filter(outliers.controls, !is.na(Specific.Group), Specific.Group == tissue)
  if (!is.na(zt)) {
    outliers.subset = filter(outliers.subset, abs(medz_value) > zt)
  }
  outliers.controls.subset = filter(outliers.controls, Gene %in% outliers.subset$Gene)
  outliers.controls.subset$Specific.Group = tissue
  return(outliers.controls.subset)
}

## Function to only keep the feature that matches the tissue group of the outliers
restrict.matching.enhancers <- function(outliers.wfeat) {
  tissue = outliers.wfeat$Specific.Group[1]
  outliers.wfeat.subset = outliers.wfeat[, c('Y','Method','Specific.Group', 'medz_value', tissue)]
  colnames(outliers.wfeat.subset)[colnames(outliers.wfeat.subset) == tissue] = 'Tissue.Enhancer'
  return(outliers.wfeat.subset)
}

################

## Read in tissue-specific outliers
data_dir = '/users/nferraro/data/goats_data/v8_data/'

fdrs = c(0.01,1e-04,1e-06)
# for (ft in fdrs) {
#   outlierfile = paste0(data_dir, 'gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.subset.Z2.FDR', ft, '.specific.RData')
#   load(outlierfile)
#   ## Input files and directories
#   featdir = paste0(baseDir, '/features_v8/byGene/500kb_genebody/')
#   outliers.controls = outliers.controls.specific
#   tissue.groups = unique(outliers.controls$Specific.Group)
#   tissue.groups = tissue.groups[!is.na(tissue.groups)]
#   
#   #### trying with all tissues together
#   outliers.wfeat = get.all.features(outliers.controls, 'MAF0-0.1', 'SNPs', featdir, '_enhancersNew_bygene.txt')
#   outliers.wfeat$feature = apply(outliers.wfeat[,9:ncol(outliers.wfeat)], 1, sum)
#   outliers.wfeat$feat.binary = (outliers.wfeat$feature > 0) + 0
#   
#   method = unique(outliers.controls$Method)
#   write.table(outliers.wfeat, file=paste0('/users/nferraro/data/goats_data/v8_data/enhancers/', method, '.outliers.Z2.FDR', ft, '.500kb.MAF0.1.wenhancers.new.V8.txt'),
#               sep='\t', quote=F, row.names=F)
#   
#   outliers.wfeat = get.all.features(outliers.controls, 'MAF0-0.1', 'indels', featdir, '_enhancersNew_bygene.txt')
#   outliers.wfeat$feature = apply(outliers.wfeat[,9:ncol(outliers.wfeat)], 1, sum)
#   outliers.wfeat$feat.binary = (outliers.wfeat$feature > 0) + 0
#   
#   write.table(outliers.wfeat, file=paste0('/users/nferraro/data/goats_data/v8_data/enhancers/', method, '.outliers.Z2.FDR', ft, '.500kb.MAF0.1.wenhancers.new.V8.indels.txt'),
#               sep='\t', quote=F, row.names=F)
#   
#   outliers.wfeat = get.all.features(outliers.controls, 'MAF0-0.1', 'HallLabSV', featdir, '_enhancersNew_bygene.txt')
#   outliers.wfeat$feature = apply(outliers.wfeat[,9:ncol(outliers.wfeat)], 1, sum)
#   outliers.wfeat$feat.binary = (outliers.wfeat$feature > 0) + 0
#   
#   write.table(outliers.wfeat, file=paste0('/users/nferraro/data/goats_data/v8_data/enhancers/', method, '.outliers.Z2.FDR', ft, '.500kb.MAF0.1.wenhancers.new.V8.HallLabSV.txt'),
#               sep='\t', quote=F, row.names=F)
# }


## PICK UP HERE
##Tissue specific analysis - START HERE

data_dir = '/users/nferraro/data/goats_data/v8_data/enhancers/'
risks = data.frame(Riskratio = numeric(), Lower = numeric(), Upper = numeric(), Pval = numeric(), MAF = numeric(), Type = character())

for (ft in fdrs) {
  knn.indel.outliers = fread(paste0(data_dir, 'KNN.outliers.Z2.FDR', ft, '.500kb.wenhancers.new.V8.indels.txt'))
  knn.snp.outliers = fread(paste0(data_dir, 'KNN.outliers.Z2.FDR', ft, '.500kb.wenhancers.new.V8.txt'))
  knn.sv.outliers = fread(paste0(data_dir, 'KNN.outliers.Z2.FDR', ft, '.500kb.wenhancers.new.V8.HallLabSV.txt'))
  knn.outliers = rbind(knn.snp.outliers, knn.indel.outliers, knn.sv.outliers)
  
  knn.specific.outliers = filter(knn.outliers, Specific.Group != '')
  knn.specific.snp.controls = filter(knn.outliers, Type == 'SNPs', is.na(Specific.Group),
                                     Gene %in% filter(knn.specific.outliers, Type == 'SNPs')$Gene)
  knn.specific.indel.controls = filter(knn.outliers, Type == 'indels', is.na(Specific.Group),
                                       Gene %in% filter(knn.specific.outliers, Type == 'indels')$Gene)
  knn.specific.sv.controls = filter(knn.outliers, Type == 'HallLabSV', Y == 0,
                                    Gene %in% filter(knn.specific.outliers, Type == 'HallLabSV')$Gene)
  knn.specific.controls = rbind(knn.specific.snp.controls, knn.specific.indel.controls, knn.specific.sv.controls) %>% mutate(Y=0)
  
  knn.specific.controls$Specific.Group = sapply(knn.specific.controls$Gene, function(x)
    knn.specific.outliers[sample(which(knn.specific.outliers$Gene == x),1),6])
  
  knn.specific = rbind(knn.specific.outliers, knn.specific.controls)
  
  tissue.feat = sapply(1:nrow(knn.specific), function(x)
    ifelse(knn.specific[x,as.character(unlist(knn.specific[x,6]))] == 1 &
             knn.specific[x,27] == 1, 1, 0))
  tissue.any = sapply(1:nrow(knn.specific), function(x)
    ifelse(knn.specific[x,as.character(unlist(knn.specific[x,6]))] == 1, 1, 0))
  knn.specific$TissFeat = tissue.feat
  knn.specific$TissAny = tissue.any
  knn.specific$TissAny = as.numeric(knn.specific$TissAny)
  knn.specific$TissAny = sapply(knn.specific$TissAny, function(x)
    ifelse(is.na(x), 0, x))
  knn.specific$TissFeat = as.numeric(knn.specific$TissFeat)
  knn.specific$TissFeat = sapply(knn.specific$TissFeat, function(x)
    ifelse(is.na(x), 0, x))
  print('Processed KNN specific')
  # Keep one line per UID when combining SNPs + indels
  knn.outliers = knn.outliers %>% group_by(Gene,Ind) %>%
    mutate(feat.sum = sum(feat.binary)) %>%
    mutate(feat.binary = ifelse(feat.sum > 0, 1, 0)) %>%
    top_n(1,Type) %>% ungroup()
  
  knn.specific = knn.specific %>% group_by(Gene,Ind) %>%
    mutate(TissFeat.sum = sum(TissFeat)) %>%
    mutate(TissAny.sum = sum(TissAny)) %>%
    mutate(TissFeat = ifelse(TissFeat.sum > 0, 1, 0)) %>%
    mutate(TissAny = ifelse(TissAny.sum > 0, 1, 0)) %>%
    top_n(1,Type) %>% ungroup()
  
  outliers = filter(knn.specific, Y == 1)
  controls = filter(knn.specific, Y == 0)
  outliers$Y = factor(outliers$Y, levels=c(0,1))
  controls$Y = factor(controls$Y, levels=c(0,1))
  counttable = rbind(table(controls$TissAny),
                     table(outliers$TissAny))
  rr = epitab(counttable,method="riskratio")$tab
  risks = rbind(risks, data.frame(Riskratio = rr[2,5], Lower = rr[2,6], Upper = rr[2,7], Pval = rr[2,8], MAF = ft, Type = 'Any'))
  counttable = rbind(table(controls$TissFeat),
                     table(outliers$TissFeat))
  print(ft)
  print(rr)
  rr = epitab(counttable,method="riskratio")$tab
  risks = rbind(risks, data.frame(Riskratio = rr[2,5], Lower = rr[2,6], Upper = rr[2,7], Pval = rr[2,8], MAF = ft, Type = 'Matched'))
}

write.table(paste0(data_dir, 'KNN.enhancer.allvars.500kb.risks.txt'),sep='\t',quote=F,row.names=F)
