library(data.table)
library(dplyr)
library(ggplot2)
library(ggthemes)
library(cowplot)

gtex_v8_figure_theme <- function() {
  return(theme(plot.title = element_text(face="plain",size=8), text = element_text(size=8),
               axis.text=element_text(size=6), panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(),panel.background = element_blank(), 
               axis.line = element_line(colour = "black"), 
               panel.border=element_blank(),
               legend.text = element_text(size=7), legend.title = element_text(size=8)))
}

data_dir = '/Users/nicoleferraro/durga_local/data/goats_data/v8_data/'
load(paste0(data_dir, 'example_genes/no_outlier_genes.RData'))
bg_genes = fread(paste0(data_dir, 'example_genes/background_genes.txt'),header=F)

te_no_outlier_genes = sapply(te_no_outlier_genes, function(x) strsplit(x, '[.]')[[1]][1])
te_no_outlier_genes = te_no_outlier_genes[which(te_no_outlier_genes %in% bg_genes$V1)]

ase_no_outlier_genes = ase_no_outlier_genes[which(ase_no_outlier_genes %in% bg_genes$V1)]

sp_no_outlier_genes = sapply(sp_no_outlier_genes, function(x) strsplit(x, '[.]')[[1]][1])
sp_no_outlier_genes = sp_no_outlier_genes[which(sp_no_outlier_genes %in% bg_genes$V1)]


medz_extreme = fread(paste0(data_dir, 'example_genes/extreme_medz_genes.txt'))
ase_extreme = fread(paste0(data_dir, 'example_genes/extreme_ase_genes.txt'))
sp_extreme = fread(paste0(data_dir, 'example_genes/extreme_splicing_genes.txt'))

medz_extreme_genes = sapply(unique(medz_extreme$Gene), function(x) strsplit(x, '[.]')[[1]][1])
sp_extreme_genes = sapply(unique(sp_extreme$CLUSTER_ID), function(x) strsplit(x, '[.]')[[1]][1])

## for entering into Panther
clip <- pipe("pbcopy", "w")                       
write.table(medz_extreme_genes, file=clip,quote=F, row.names=F, col.names=F)                               
close(clip)

## plot GO enrichments for no outlier genes
gcols = c('GO_term', 'BG', 'Test_genes', 'BG_expected', 'Test_over_under', 'Fold_Enrichment', 'Pval', 'FDR')
te_no_outlier_GO = fread(paste0(data_dir, 'example_genes/te_no_outlier_genes_panther_GO_bp.txt'), sep='\t', fill=T, col.names=gcols)
ase_no_outlier_GO = fread(paste0(data_dir, 'example_genes/ase_no_outlier_genes_panther_GO_slim_bp.txt'), sep='\t', fill=T, col.names=gcols)
sp_no_outlier_GO = fread(paste0(data_dir, 'example_genes/splicing_no_outlier_genes_panther_GO_slim_bp.txt'), sep='\t', fill=T, col.names=gcols)

ccols = c('grey', alpha('maroon4',0.7))
names(ccols) = c('no', 'yes')

te_GO_sum = te_no_outlier_GO %>% top_n(10,-FDR) %>% arrange(by=-FDR)
te_GO_sum$GO_term = sapply(te_GO_sum$GO_term, function(x) strsplit(x, ' [(]GO')[[1]][1])
te_GO_sum$GO_term = factor(te_GO_sum$GO_term, levels=unique(te_GO_sum$GO_term))
te_GO_sum$Significant = sapply(te_GO_sum$FDR, function(x) ifelse(x < 0.05, 'yes', 'no'))

te_plot = ggplot(te_GO_sum, aes(x=GO_term,y=-log10(FDR))) + 
  geom_bar(stat='identity',aes(fill=Significant),width=0.5) + 
  theme_bw() + coord_flip() + xlab('') + 
  scale_fill_manual(values=ccols) + guides(fill=F) +
  ggtitle('eOutlier (n=261)') + gtex_v8_figure_theme()

ase_GO_sum = ase_no_outlier_GO %>% top_n(10,-FDR) %>% arrange(by=-FDR)
ase_GO_sum$GO_term = sapply(ase_GO_sum$GO_term, function(x) strsplit(x, ' [(]GO')[[1]][1])
ase_GO_sum$GO_term = factor(ase_GO_sum$GO_term, levels=unique(ase_GO_sum$GO_term))
ase_GO_sum$Significant = sapply(ase_GO_sum$FDR, function(x) ifelse(x < 0.05, 'yes', 'no'))

ase_plot = ggplot(ase_GO_sum, aes(x=GO_term,y=-log10(FDR))) + 
  geom_bar(stat='identity', aes(fill=Significant),width=0.5) + theme_bw() + coord_flip() + xlab('') + 
  scale_fill_manual(values=ccols) + guides(fill=F) +
  ggtitle('aseOutlier (n=11573)') + gtex_v8_figure_theme()

sp_GO_sum = sp_no_outlier_GO %>% top_n(10,-FDR) %>% arrange(by=-FDR)
sp_GO_sum$GO_term = sapply(sp_GO_sum$GO_term, function(x) strsplit(x, ' [(]GO')[[1]][1])
sp_GO_sum$GO_term = factor(sp_GO_sum$GO_term, levels=unique(sp_GO_sum$GO_term))
sp_GO_sum$Significant = sapply(sp_GO_sum$FDR, function(x) ifelse(x < 0.05, 'yes', 'no'))

sp_plot = ggplot(sp_GO_sum, aes(x=GO_term,y=-log10(FDR))) + 
  geom_bar(stat='identity', aes(fill=Significant),width=0.5) + theme_bw() + coord_flip() + xlab('') + 
  scale_fill_manual(values=ccols) + guides(fill=F) +
  ggtitle('sOutlier (n=9548)') + gtex_v8_figure_theme()

first_row = plot_grid(te_plot, ase_plot, sp_plot, nrow = 1, axis='tb')

te_ex_outlier_GO = fread(paste0(data_dir, 'example_genes/te_extreme_outlier_genes_panther_GO_bp.txt'), sep='\t', fill=T, col.names=gcols)
ase_ex_outlier_GO = fread(paste0(data_dir, 'example_genes/ase_extreme_outlier_genes_panther_GO_bp.txt'), sep='\t', fill=T, col.names=gcols)
sp_ex_outlier_GO = fread(paste0(data_dir, 'example_genes/splicing_extreme_outlier_genes_panther_GO_bp.txt'), sep='\t', fill=T, col.names=gcols)

te_GO_sum = te_ex_outlier_GO %>% top_n(10,-Pval) %>% arrange(by=-Pval)
te_GO_sum$GO_term = sapply(te_GO_sum$GO_term, function(x) strsplit(x, ' [(]GO')[[1]][1])
te_GO_sum$GO_term = factor(te_GO_sum$GO_term, levels=unique(te_GO_sum$GO_term))
te_GO_sum$Significant = sapply(te_GO_sum$FDR, function(x) ifelse(x < 0.05, 'yes', 'no'))

te_plot2 = ggplot(te_GO_sum, aes(x=GO_term,y=-log10(Pval))) + 
  geom_bar(stat='identity',aes(fill=Significant),width=0.5) + theme_bw() + coord_flip() + xlab('') + 
  scale_fill_manual(values=ccols) + guides(fill=F) +
  ggtitle('eOutlier (n=127)') + gtex_v8_figure_theme()

ase_GO_sum = ase_ex_outlier_GO %>% top_n(10,-FDR) %>% arrange(by=-FDR)
ase_GO_sum$GO_term = sapply(ase_GO_sum$GO_term, function(x) strsplit(x, ' [(]GO')[[1]][1])
ase_GO_sum$GO_term = factor(ase_GO_sum$GO_term, levels=unique(ase_GO_sum$GO_term))
ase_GO_sum$Significant = sapply(ase_GO_sum$FDR, function(x) ifelse(x < 0.05, 'yes', 'no'))

ase_plot2 = ggplot(ase_GO_sum, aes(x=GO_term,y=-log10(FDR))) + 
  geom_bar(stat='identity',aes(fill=Significant),width=0.5) + theme_bw() + coord_flip() + xlab('') + 
  scale_fill_manual(values=ccols) + guides(fill=F) +
  ggtitle('aseOutlier (n=261)') + gtex_v8_figure_theme()

sp_GO_sum = sp_ex_outlier_GO %>% top_n(10,-FDR) %>% arrange(by=-FDR)
sp_GO_sum$GO_term = sapply(sp_GO_sum$GO_term, function(x) strsplit(x, ' [(]GO')[[1]][1])
sp_GO_sum$GO_term = factor(sp_GO_sum$GO_term, levels=unique(sp_GO_sum$GO_term))
sp_GO_sum$Significant = sapply(sp_GO_sum$FDR, function(x) ifelse(x < 0.05, 'yes', 'no'))

sp_plot2 = ggplot(sp_GO_sum, aes(x=GO_term,y=-log10(FDR))) + 
  geom_bar(stat='identity',aes(fill=Significant),width=0.5) + theme_bw() + coord_flip() + xlab('') + 
  scale_fill_manual(values=ccols) + guides(fill=F) +
  ggtitle('sOutlier (n=389)') + gtex_v8_figure_theme()

sfig = plot_grid(te_plot, ase_plot, sp_plot, te_plot2, ase_plot2, sp_plot2, nrow = 2, axis='tb', align='hv')
ggsave(sfig, file=paste0(fig_dir, 'sfig_go_enrich.pdf'), width=7.2, height=4, units="in", scale=2)




