#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
require(dplyr)
library(reshape)
library(scales)
require(RColorBrewer)

data_dir = '/Users/nicoleferraro/durga_local/data/goats_data/v8_data/'

medz_annots = fread(paste0(data_dir, 'gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.medz.variantAnnotations.txt'),data.table=F)
medz_annots = filter(medz_annots,sv_v7 == 1)

splice_annot = fread(paste0(data_dir, 'splicing/gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.splice.variantAnnotations.filtered.txt'),data.table=F)
splice_annot = filter(splice_annot,sv_v7==1)

ase_annot = fread(paste0(data_dir,'ASE/gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.ase.variantAnnotations.txt'),data.table=F)
ase_annot = filter(ase_annot,sv_v7==1)

wgs_inds = fread('/Users/nicoleferraro/durga_local/data/goats_data/v8_data/gtexv8_wgs_sv_inds.txt',header=F)$V1
wgs_inds = sapply(wgs_inds, function(x) strsplit(x,'GTEX-')[[1]][2])

medz_annots$indiv_id = factor(medz_annots$indiv_id,levels=wgs_inds)
splice_annot$indiv_id = factor(splice_annot$indiv_id,levels=wgs_inds)
ase_annot$indiv_id = factor(ase_annot$indiv_id,levels=wgs_inds)

num_medz_per_ind = as.data.frame(table(filter(medz_annots,abs(MedZ)>3)$indiv_id)) %>%
  mutate(Category='Outliers') %>% mutate(Type = 'TotalExp')
num_splice_per_ind = as.data.frame(table(filter(splice_annot,value < 0.0027)$indiv_id)) %>%
  mutate(Category='Outliers') %>% mutate(Type = 'Splicing')
num_ase_per_ind = as.data.frame(table(filter(ase_annot,DOT.score < 0.0027)$indiv_id))  %>%
  mutate(Category='Outliers') %>% mutate(Type = 'ASE')

num_medz_RVs_per_ind = as.data.frame(table(filter(medz_annots,abs(MedZ)>3,variant_cat != 'no_variant')$indiv_id)) %>%
  mutate(Category='Outliers with RV') %>% mutate(Type = 'TotalExp')
num_splice_RVs_per_ind = as.data.frame(table(filter(splice_annot,value < 0.0027,variant_cat != 'no_variant')$indiv_id)) %>%
  mutate(Category='Outliers with RV') %>% mutate(Type = 'Splicing')
num_ase_RVs_per_ind = as.data.frame(table(filter(ase_annot,DOT.score < 0.0027,variant_cat != 'no_variant')$indiv_id)) %>%
  mutate(Category='Outliers with RV') %>% mutate(Type = 'ASE')

num_overall = rbind(num_medz_per_ind,num_splice_per_ind,num_ase_per_ind,
                    num_medz_RVs_per_ind,num_splice_RVs_per_ind, num_ase_RVs_per_ind)

tcols = brewer.pal(9,'BuPu')[c(3,5,7)]
ggplot(num_overall, aes(x=Category,y=Freq,Group=Type)) +
  geom_boxplot(aes(fill=Type)) + theme_bw() +
  scale_fill_manual(values=tcols) +
  xlab('') + ylab('Number per ind') +
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=20),
        legend.title = element_blank(),
        legend.text = element_text(size=18)) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

wdata = fread(paste0(data_dir, 'coloc/watershed_posteriors.txt'),data.table=F)










