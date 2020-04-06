## check variants in TOPMed

library(ggplot2)
library(data.table)
library(dplyr)
library(RColorBrewer)
library(cowplot)


gtex_v8_figure_theme <- function() {
  return(theme(plot.title = element_text(face="plain",size=8), text = element_text(size=8),
               axis.text=element_text(size=7), panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(),panel.background = element_blank(), 
               panel.border = element_blank(),
               axis.line = element_line(colour = "black"), 
               legend.text = element_text(size=7), legend.title = element_text(size=8)))
}

data_dir = '/Users/nicoleferraro/durga_local/data/goats_data/v8_data/enrichments/windows/'

exp_counts = fread(paste0(data_dir, 'S07_Check_Spatialv2/Expression_counts.txt')) %>% mutate(Type = 'eOutliers')
ase_counts = fread(paste0(data_dir, 'v8_vg/S07_Check_Spatialv2/ASE_v8_counts.txt')) %>% mutate(Type = 'aseOutliers')
sp_counts = fread(paste0(data_dir, 'S07_Check_Spatialv2/Splicing_CoClusterFiltered_counts.txt')) %>% mutate(Type = 'sOutliers')
all_counts = rbind(exp_counts, ase_counts, sp_counts)
colnames(all_counts) = c('Window', 'Observed', 'Expected', 'Type')
all_counts_melted = melt(all_counts, id.vars=c('Window', 'Type'))
all_counts_melted$Type = factor(all_counts_melted$Type, levels=c('eOutliers', 'aseOutliers', 'sOutliers'))
sfig_A = ggplot(all_counts_melted, aes(x=Window, y=value, Group=variable)) +
  geom_bar(stat='identity',aes(fill=variable),position='dodge') + 
  gtex_v8_figure_theme() + facet_wrap(.~Type) + xlab('Window size: log10(bp)') +
  ylab('Outlier pairs within window') +
  theme(legend.title=element_blank(),
        strip.background = element_blank(),
        legend.position=c(0.1,0.8),
        legend.key.height = unit(0.05, "cm"),
        panel.spacing = unit(2, "lines"))

exp_enrich = fread(paste0(data_dir, 'S07_Check_Spatialv2/Expression_Enrichment.txt')) %>% mutate(Type = 'eOutliers')
ase_enrich = fread(paste0(data_dir, 'v8_vg/S07_Check_Spatialv2/ASE_v8_Enrichment.txt')) %>% mutate(Type = 'aseOutliers')
sp_enrich = fread(paste0(data_dir, 'S07_Check_Spatialv2/Splicing_CoClusterFiltered_Enrichment.txt')) %>% mutate(Type = 'sOutliers')
all_enrich = rbind(exp_enrich, ase_enrich, sp_enrich)
colnames(all_enrich) = c('Enrichment', 'Type')
all_enrich$Window = rep(exp_counts$x_log10bp_WindowSize, 3)
all_enrich$LogRatio = log2(all_enrich$Enrichment)
all_enrich$LogRatio = sapply(all_enrich$LogRatio, function(x) ifelse(is.infinite(x), 0, x))
all_enrich$Type = factor(all_enrich$Type, levels=c('eOutliers', 'aseOutliers', 'sOutliers'))
sfig_B = ggplot(all_enrich, aes(x=Window, y=LogRatio)) +
  geom_bar(stat='identity') + geom_hline(yintercept=0, color='lightgrey') +
  gtex_v8_figure_theme() + facet_wrap(.~Type) + xlab('Window size: log10(bp)') +
  ylab('Expected / Observed') +
  theme(strip.background = element_blank(),
        panel.spacing = unit(2, "lines"))

sp_enrich = fread(paste0(data_dir, 'S07_Check_Spatialv2/Splicing_Enrichment.txt')) %>% mutate(Type = 'sOutliers')
colnames(sp_enrich) = c('Enrichment', 'Type')
sp_enrich$Window = exp_counts$x_log10bp_WindowSize
sp_enrich$LogRatio = log2(sp_enrich$Enrichment)
sp_enrich$LogRatio = sapply(sp_enrich$LogRatio, function(x) ifelse(is.infinite(x), 0, x))
sfig_C = ggplot(sp_enrich, aes(x=Window, y=LogRatio)) +
  geom_bar(stat='identity') + geom_hline(yintercept=0, color='lightgrey') +
  ggtitle('sOutlier enrichments w/o filtering') +
  gtex_v8_figure_theme() + xlab('Window size: log10(bp)') +
  ylab('Expected / Observed') +
  theme(strip.background = element_blank())

all_files = dir(data_dir, 'ase.te.pair.enrichments.*.txt', full=T)
all.rrs = do.call(rbind,lapply(all_files, function(x) fread(x)))
all.rrs = filter(all.rrs, !(VCat %in% c(13,14))) %>% filter(!is.na(Upper)) %>% filter(!is.na(Lower))
all.rrs$Window = factor(all.rrs$Window, levels=c('100kbp', '500kbp', '1000kbp', '5000kbp'))
sfig_D = ggplot(all.rrs %>% filter(Window == '100kbp', DataType == 'TE', Cat %in% c('CNV', 'DUP', 'TSS')), aes(x=Cat, y=Riskratio)) + 
  geom_point(size=3) + xlab('') + ylab('Relative risk') +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width=0) +
  gtex_v8_figure_theme() + scale_y_log10() + 
  ggtitle('eOutlier pairs') +
  geom_hline(yintercept=1, linetype="dashed", colour="darkgrey")

third_row <- plot_grid(sfig_C, sfig_D, labels=c('C', 'D'), ncol=2, align='hv', axis='tblr')
sfig <- plot_grid(sfig_A, sfig_B, third_row, labels=c('A', 'B', NULL, NULL), nrow=3)

fig_dir = '/Users/nicoleferraro/durga_local/data/goats_data/v8_data/paper_figures/'
ggsave(sfig, file=paste0(fig_dir, 'sfig_window_v4.pdf'), width=7.2, height=9, units="in")



