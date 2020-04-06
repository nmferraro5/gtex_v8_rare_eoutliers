#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
require(dplyr)
library(reshape)
library(scales)
require(RColorBrewer)

data_dir = '/Users/nicoleferraro/durga_local/data/goats_data/v8_data/'
data_dir = '/users/nferraro/data/goats_data/v8_data/'

pvals = c(0.05, 0.01, 0.001, 0.0001, 0.00001, 0.0000001)
zvals = sapply(pvals, function(x) qnorm(1 - x))
get_bin <- function(mval) {
  if (abs(mval) > zvals[6]) {
    return('0~1e-07')
  } else if (abs(mval) > zvals[5]) {
    return('1e-07~1e-05')
  } else if (abs(mval) > zvals[4]) {
    return('1e-05~1e-04')
  } else if (abs(mval) > zvals[3]) {
    return('1e-04~1e-03')
  } else if (abs(mval) > zvals[2]) {
    return('1e-03~1e-02')
  } else if (abs(mval) > zvals[1]) {
    return('1e-02~5e-02')
  } else {
    return('nonOutlier')
  }
}

get_pbin <- function(mval) {
  if (mval < pvals[6]) {
    return('0~1e-07')
  } else if (mval < pvals[5]) {
    return('1e-07~1e-05')
  } else if (mval < pvals[4]) {
    return('1e-05~1e-04')
  } else if (mval < pvals[3]) {
    return('1e-04~1e-03')
  } else if (mval < pvals[2]) {
    return('1e-03~1e-02')
  } else if (mval < pvals[1]) {
    return('1e-02~5e-02')
  } else {
    return('nonOutlier')
  }
}

get_fc_bin <- function(afc) {
  if (is.na(afc)) {
    return('nonOutlier')
  }
  else if (abs(afc) > 2) {
    return('2~Inf')
  } else if (abs(afc) > 1.5) {
    return('1.5~2')
  } else if (abs(afc) > 1) {
    return('1~1.5')
  } else if (abs(afc) > 0.5) {
    return('0.5~1')
  } else {
    return('0~0.5')
  }
}

var_data = fread(paste0(data_dir, 'gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.medz.variantAnnotations.filtered.txt'),data.table=F)
var_data = var_data %>% select(indiv_id,gene_id,MedZ,Y,indexInRow,groupCount,priority,tier2,af_gtex,categoryOutlier,variant_cat,LoF,sv_v7,af_gnomad,af_bin,medz_bin,medz_bin2,color_R,color_G,color_B)
var_data$medz_bin3 = sapply(var_data$MedZ, function(x) get_bin(x))
all_var_data = var_data
var_data = filter(all_var_data, medz_bin3 == 'nonOutlier' | MedZ > 0)

ase_logfc = fread(paste0(data_dir, 'ASE/MEDIAN_dots_and_aFC.tsv'),data.table=F)
ase_logfc$indiv_id = sapply(ase_logfc$SampleName, function(x) strsplit(x, 'GTEX-')[[1]][2])
colnames(ase_logfc)[1] = 'gene_id'
var_data = merge(var_data, ase_logfc %>% dplyr::select(gene_id,indiv_id,median_aFC), by=c('gene_id','indiv_id'),all.x=T)
var_data$FC_bin = sapply(var_data$median_aFC, function(x) get_fc_bin(x))

var_data = fread(paste0(data_dir, 'ASE/gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.ase.variantAnnotations.v8.vg.txt'))
var_data = var_data %>% select(indiv_id,gene_id,value,pval_bin,tier2,variant_cat,sv_v7,af_gnomad,variant_color1,variant_color2,variant_color3) # previously color_R, color_G, color_B
var_data = filter(var_data,sv_v7==1)
new_cats = sapply(1:nrow(var_data), function(x) ifelse(var_data$variant_cat[x] == 'splice', var_data$tier2[x], var_data$variant_cat[x]))
var_data$variant_cat = new_cats
var_data = filter(var_data, !is.na(value))
pbins = sapply(var_data$value, function(x) get_pbin(x))
var_data$pval_bin = pbins
plot_medz_bin = var_data %>% 
  group_by(pval_bin) %>%
  mutate(NumBin = n()) %>% ungroup() %>%
  group_by(pval_bin,variant_cat) %>%
  mutate(NumCat = n()/NumBin) %>% sample_n(1) %>% ungroup()

#plot_medz_bin = plot_medz_bin %>% mutate(CatCol = rgb(color_R, color_G, color_B))
plot_medz_bin = plot_medz_bin %>% mutate(CatCol = rgb(variant_color1, variant_color2, variant_color3))

plot_medz_bin = var_data %>%
 group_by(medz_bin3) %>%
 mutate(NumBin = n()) %>% ungroup() %>%
 group_by(medz_bin3,variant_cat) %>%
 mutate(NumCat = n()/NumBin) %>% sample_n(1) %>% ungroup()

plot_medz_bin = plot_medz_bin %>% filter(!(medz_bin %in% c(NA,"0.000~0.002"))) %>%
 mutate(CatCol = rgb(color_R, color_G, color_B))

#plot_cols = unique(plot_medz_bin$CatCol)
#names(plot_cols) = unique(plot_medz_bin$variant_cat)
print(unique(plot_medz_bin$variant_cat))
#plot_cols = c(unique(plot_medz_bin$CatCol)[1:10],'#3e6690','#33669a',unique(plot_medz_bin$CatCol)[11:16]) # ends at 9 for splice
plot_cols = c(unique(plot_medz_bin$CatCol)[1:11],'#3e6690','#33669a',unique(plot_medz_bin$CatCol)[12:16]) # ends at 9 for splice
names(plot_cols) = unique(plot_medz_bin$variant_cat)

plot_medz_bin$pval_bin = factor(plot_medz_bin$pval_bin,
                                levels=c('0~1e-07','1e-07~1e-05','1e-05~1e-04','1e-04~1e-03','1e-03~1e-02','1e-02~5e-02','nonOutlier'))
plot_medz_bin$medz_bin3 = factor(plot_medz_bin$medz_bin3,
                                levels=c('0~1e-07','1e-07~1e-05','1e-05~1e-04','1e-04~1e-03','1e-03~1e-02','1e-02~5e-02','nonOutlier'))

plot_medz_bin$FC_bin = factor(plot_medz_bin$FC_bin, 
                              levels=c('2~Inf', '1.5~2', '1~1.5', '0.5~1', '0~0.5', 'nonOutlier'))
#plot_medz_bin$medz_bin = factor(plot_medz_bin$medz_bin,
#                                levels=c('-Inf~-5','-5~-4','-4~-3','-3~-2','-2~-1','-1~1','1~2','2~3','3~4','4~5','5~Inf'))
#plot_medz_bin$variant_cat = factor(plot_medz_bin$variant_cat,
#                                   levels=c('no_variant','other_noncoding','coding','conserved_noncoding','TSS','stop','frameshift','splice','TE','INV','BND','DEL','CNV','DUP'))
plot_medz_bin$variant_cat = factor(plot_medz_bin$variant_cat,
                                   levels=c('no_variant','other_noncoding','coding','conserved_noncoding','TSS','stop','frameshift','splice_region_variant','splice_donor_variant','splice_acceptor_variant','TE','INV','BND','DEL','CNV','DUP'))

vplot = ggplot(plot_medz_bin, aes(x=pval_bin,y=NumCat,Group=variant_cat)) +
  geom_bar(aes(fill=variant_cat),color='black',stat='identity') + 
  theme_bw() + coord_flip() + ylab('') + xlab('') +
  scale_fill_manual(values=plot_cols) + 
  ggtitle('Total Expression - Over') + xlab('P-value bin') + 
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.text=element_text(size=16),
        title=element_text(size=22),
        legend.title=element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

pdf('/Users/nicoleferraro/durga_local/data/goats_figures/all.variant.barplots.pdf', height=13,width=8)
multiplot(te_over_plot, te_under_plot, splice_plot, ase_plot, cols=1)
dev.off()

save(plot_ase_bin,vplot,file=paste0(data_dir,'ase_v8_vg_variant_annotation_plot.RData'))
     
     
