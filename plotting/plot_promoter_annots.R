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

get_pbin <- function(pval) {
  if (is.na(pval)) {
    return('nonOutlier')
  }
  else if (pval < pvals[6]) {
    return('0~1e-07')
  } else if (pval < pvals[5]) {
    return('1e-07~1e-05')
  } else if (pval < pvals[4]) {
    return('1e-05~1e-04')
  } else if (pval < pvals[3]) {
    return('1e-04~1e-03')
  } else if (pval < pvals[2]) {
    return('1e-03~1e-02')
  } else if (pval < pvals[1]) {
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

medz_data = fread(paste0(data_dir, 'gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.medz.variantAnnotations.withPromoter.filtered.outlierGenes.txt'),data.table=F)
ase_data = fread(paste0(data_dir, 'ASE/gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.ase.variantAnnotations.withPromoter.filtered.txt'))

medz_tss_data = filter(medz_data, variant_cat == 'TSS')
#medz_tss_data$medz_bin3 = sapply(medz_tss_data$MedZ, function(x) get_bin(x))
medz_tss_data$medz_bin = sapply(medz_tss_data$MedZ, function(x) ifelse(x < -3, 'Under', 
                                                                      ifelse(x < 3, 'Control','Over')))

ase_data$pval_bin = sapply(ase_data$value, function(x) get_pbin(x))
ase_logfc = fread(paste0(data_dir, 'ASE/MEDIAN_dots_and_aFC.tsv'),data.table=F)
ase_logfc$indiv_id = sapply(ase_logfc$SampleName, function(x) strsplit(x, 'GTEX-')[[1]][2])
colnames(ase_logfc)[1] = 'gene_id'
ase_data = merge(ase_data, ase_logfc %>% dplyr::select(gene_id,indiv_id,median_aFC), by=c('gene_id','indiv_id'),all.x=T)
ase_data$FC_bin = sapply(ase_data$median_aFC, function(x) get_fc_bin(x))
ase_tss_data = filter(ase_data, variant_cat == 'TSS')

plot_medz_bin = medz_tss_data %>% 
  group_by(medz_bin) %>%
  mutate(NumBin = n()) %>% ungroup() %>%
  group_by(medz_bin,promoter_motif) %>%
  mutate(NumCat = n()/NumBin) %>% sample_n(1) %>% ungroup()

plot_medz_bin = plot_medz_bin %>% mutate(CatCol = rgb(variant_color1, variant_color2, variant_color3))

plot_ase_bin = ase_tss_data %>% 
  group_by(FC_bin) %>%
  mutate(NumBin = n()) %>% ungroup() %>%
  group_by(FC_bin,promoter_motif) %>%
  mutate(NumCat = n()/NumBin) %>% sample_n(1) %>% ungroup()

plot_ase_bin = plot_ase_bin %>% mutate(CatCol = rgb(variant_color1, variant_color2, variant_color3))


plot_ase_bin$pval_bin = factor(plot_ase_bin$pval_bin,
                                levels=c('0~1e-07','1e-07~1e-05','1e-05~1e-04','1e-04~1e-03','1e-03~1e-02','1e-02~5e-02','nonOutlier'))
plot_medz_bin$medz_bin3 = factor(plot_medz_bin$medz_bin3,
                                levels=c('0~1e-07','1e-07~1e-05','1e-05~1e-04','1e-04~1e-03','1e-03~1e-02','1e-02~5e-02','nonOutlier'))

plot_medz_bin$promoter_motif = factor(plot_medz_bin$promoter_motif, 
                                      levels=c('no_motif', 'other_motif', 'Cmyc', 'CTCF', 'E2F4', 'E2F6', 'ELF1', 'Gabp', 'Nrf1', 'Nrsf', 'PU1', 'SP1', 'Srf', 'Tr4', 'USF1', 'Yy1'))
pcols = c('white', 'lightgrey', brewer.pal(7,'BrBG'),brewer.pal(7,'RdYlBu'))
plot_medz_bin$medz_bin = factor(plot_medz_bin$medz_bin, levels=c('Under', 'Control', 'Over'))
names(pcols) = levels(plot_medz_bin$promoter_motif)

plot_ase_bin$promoter_motif = factor(plot_ase_bin$promoter_motif, 
                                      levels=c('no_motif', 'other_motif', 'Cmyc', 'CTCF', 'E2F4', 'E2F6', 'ELF1', 'Gabp', 'Nrf1', 'Nrsf', 'PU1', 'SP1', 'Srf', 'Tr4', 'USF1', 'Yy1'))
pcols = c('white', 'lightgrey', brewer.pal(7,'BrBG'),brewer.pal(7,'RdYlBu'))
plot_ase_bin$FC_bin = factor(plot_ase_bin$FC_bin, levels=c('2~Inf', '1.5~2', '1~1.5', '0.5~1', '0~0.5', 'nonOutlier'))


vplot = ggplot(plot_medz_bin, aes(x=medz_bin,y=NumCat,Group=promoter_motif)) +
  geom_bar(aes(fill=promoter_motif),color='black',stat='identity') + 
  theme_bw() + ylab('') + xlab('') +
  ylim(c(0,0.4)) + guides(fill=F) +
  scale_fill_manual(values=pcols) +
  ylab('Proportion rare variants in promoter') +
  ggtitle('Total Expression') + xlab('') + 
  theme(axis.text=element_text(size=20),
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

save(plot_medz_bin,vplot,file=paste0(data_dir,'medz_variant_annotation_plot.RData'))
