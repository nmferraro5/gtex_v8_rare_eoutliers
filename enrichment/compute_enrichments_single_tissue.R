#!/usr/bin/env Rscript

rm(list = ls())

library(data.table)
library(dplyr)
library(stringr)
library(epitools)

## Script for combining outlier calls and features then computing enrichments.

baseDir = Sys.getenv('RAREDIR')

data_dir = '/users/nferraro/data/goats_data/v8_data/'
exp_single_files = dir(paste0(baseDir, '/preprocessing_v8/PEER_v8'), '*.peer.v8ciseQTLs.ztrans.txt', full=T)
print(exp_single_files[1])
rare_snps = fread(paste0(data_dir, 'all_genes_with_rare_SNPs.txt'), header=F)
rare_indels = fread(paste0(data_dir, 'all_genes_with_rare_indels.txt'), header=F)
rare_svs = fread(paste0(data_dir, 'all_genes_with_rare_HallLabSV.txt'), header=F)
colnames(rare_snps) = c('Ind','Gene')
colnames(rare_indels) = c('Ind','Gene')
colnames(rare_svs) = c('Ind','Gene')
#rare_snps$Gene = sapply(rare_snps$Gene, function(x) strsplit(x, '[.]')[[1]][1])
#rare_indels$Gene = sapply(rare_indels$Gene, function(x) strsplit(x, '[.]')[[1]][1])
#rare_svs$Gene = sapply(rare_svs$Gene, function(x) strsplit(x, '[.]')[[1]][1])
rare_snps = rare_snps %>% mutate(UID = paste(Gene,Ind,sep='_')) %>% select(UID)
rare_indels = rare_indels %>% mutate(UID = paste(Gene,Ind,sep='_')) %>% select(UID)
rare_svs = rare_svs %>% mutate(UID = paste(Gene,Ind,sep='_'))

#pvals = c(0.05, 0.01, 0.001, 0.0001, 0.00001, 0.0000001, 0.000000001, 0.00000000001, 0.0000000000001)
zvals = c(1,2,3,4,5,6,7,8,9,10)
pvals = sapply(zvals, function(x) 2*pnorm(-abs(x)))

#ase_data = fread(paste0(baseDir, '/data_v8/ase/combined.uncorrected.unfiltered.ad.scores.all.txt')) %>% mutate(UID = paste(V1,variable,sep='_'))
#colnames(ase_data)[1:2] = c('Gene','Ind')
#ase_inds = fread(paste0(data_dir, 'ASE/ase_inds_update.txt'))
#ase_data = filter(ase_data, Ind %in% rare_svs$Ind)
#splice_data = fread(paste0(baseDir, '/data_v8/splicing/all_tissue_covariate_method_none_no_global_outliers_ea_only_emperical_pvalue_gene_level.txt'))
#splice_melted = melt(splice_data)
#colnames(splice_melted)[3] = 'Ind'
#splice_melted = splice_melted %>% filter(Ind %in% rare_svs$Ind)
#print(splice_melted[1,])
goutliers = fread(paste0(data_dir, 'gtexV8_global_outliers_medz3_iqr.txt'),header=F)
estims = data.frame(Riskratio = numeric(), Lower = numeric(), Upper = numeric(), Pval = numeric(), ZT = numeric(), Type = character(), Tissue = character())
for (zt in zvals) {
    print(zt)
    for (ef in exp_single_files) {
        gdata = melt(fread(ef))
        tissue = strsplit(ef, '/')[[1]][8]
        tissue = strsplit(tissue, '.peer')[[1]][1]
        print(tissue)
        colnames(gdata) = c('Gene','Ind','Zscore')
        #gdata$Pval = 2*pnorm(-abs(gdata$Zscore))
        gdata = gdata %>% filter(!(Ind %in% goutliers$V1))
        #outliers = filter(gdata, Pval < pt) %>% select(-Zscore,-Pval) %>% mutate(UID = paste(Gene,Ind,sep='_'))
        #controls = filter(gdata, Pval > pt, Gene %in% outliers$Gene) %>% select(-Zscore,-Pval) %>% mutate(UID = paste(Gene,Ind,sep='_')) %>% select(UID)
        outliers = filter(gdata, abs(Zscore) > zt) %>% select(-Zscore) %>% mutate(UID = paste(Gene,Ind,sep='_'))
        controls = filter(gdata, abs(Zscore) < zt, Gene %in% outliers$Gene) %>% select(-Zscore) %>% mutate(UID = paste(Gene,Ind,sep='_')) %>% select(UID)
        yy_snps = length(which(outliers$UID %in% rare_snps$UID))
        ny_snps = length(which(controls$UID %in% rare_snps$UID))
        nn_snps = nrow(controls) - ny_snps
        yn_snps = nrow(outliers) - yy_snps
        yy_indels = length(which(outliers$UID %in% rare_indels$UID))
        ny_indels = length(which(controls$UID %in% rare_indels$UID))
        nn_indels = nrow(controls) - ny_indels
        yn_indels = nrow(outliers) - yy_indels
        yy_svs = length(which(outliers$UID %in% rare_svs$UID))
        ny_svs = length(which(controls$UID %in% rare_svs$UID))
        nn_svs = nrow(controls) - ny_svs
        yn_svs = nrow(outliers) - yy_svs
        counttable = rbind(c(nn_snps,ny_snps),c(yn_snps,yy_snps))
        tryCatch({
            rr = epitab(counttable,method="riskratio")$tab
            estims = rbind(estims, data.frame(Riskratio = rr[2,5],
                                                Lower = rr[2,6],
                                                Upper = rr[2,7],
                                                Pval = rr[2, 8],
                                                ZT = zt,
                                                Type = 'SNPs',
                                                NOutliers = nrow(outliers),
                                                Tissue = tissue))
        }, error = function(e) { print(counttable) })
        counttable = rbind(c(nn_indels,ny_indels),c(yn_indels,yy_indels))
        tryCatch({
            rr = epitab(counttable,method="riskratio")$tab
            estims = rbind(estims, data.frame(Riskratio = rr[2,5],
                                                Lower = rr[2,6],
                                                Upper = rr[2,7],
                                                Pval = rr[2, 8],
                                                ZT = zt,
                                                Type = 'indels',
                                                NOutliers = nrow(outliers),
                                                Tissue = tissue))
        }, error = function(e) { print(counttable) })
        counttable = rbind(c(nn_svs,ny_svs),c(yn_svs,yy_svs))
        tryCatch({
            rr = epitab(counttable,method="riskratio")$tab
            estims = rbind(estims, data.frame(Riskratio = rr[2,5],
                                                Lower = rr[2,6],
                                                Upper = rr[2,7],
                                                Pval = rr[2, 8],
                                                ZT = zt,
                                                Type = 'SVs',
                                                NOutliers = nrow(outliers),
                                                Tissue = tissue))
        }, error = function(e) { print(counttable) })
    }
}

write.table(estims, file=paste0(data_dir, 'enrichments/all.single.tissue.expression.relative.risks.txt'), sep='\t', quote=F, row.names=F)
