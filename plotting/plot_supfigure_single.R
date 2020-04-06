### analysis for supplemental figure 2 in section 1 for TE - comparison of outlier types

library(gplots)
library(ggplot2)
library(ggthemes)
library(ggridges)
library(data.table)
library(dplyr)
library(RColorBrewer)
library(cowplot)

gtex_v8_figure_theme <- function() {
  return(theme(plot.title = element_text(face="plain",size=8), text = element_text(size=8),
               axis.text=element_text(size=7), panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(),panel.background = element_blank(), 
               axis.line = element_line(colour = "black"), 
               panel.border=element_blank(),
               legend.text = element_text(size=7), legend.title = element_text(size=8)))
}

data_dir = '/Users/nicoleferraro/durga_local/data/goats_data/v8_data/'

single_enrich = fread(paste0(data_dir, 'enrichments/all.single.tissue.expression.relative.risks.txt'))
single_enrich$Type = sapply(single_enrich$Type, function(x) ifelse(x == 'SNPs', 'SNVs', x))
single_enrich$ZT = factor(single_enrich$ZT, levels=unique(single_enrich$ZT))
single_enrich$Type = factor(single_enrich$Type, levels=c('SNVs', 'indels', 'SVs'))
senrich = ggplot(single_enrich, aes(x=ZT, y=Riskratio)) + geom_boxplot() + facet_grid(Type~., scales='free') +
  geom_hline(yintercept=1,color='firebrick4') + theme_bw() +
  scale_y_log10() + xlab('|Z-score| threshold') + ylab('Relative risk') +
  gtex_v8_figure_theme()

sfig_A = ggplot(single_enrich %>% filter(Type == 'SNVs'), aes(x=ZT, y=Riskratio)) + geom_boxplot() + 
  xlab('|Z-score| threshold') + ylab('Relative risk') + ggtitle('SNVs') +
  gtex_v8_figure_theme() + geom_hline(yintercept=1, color='darkgrey')

sfig_B = ggplot(single_enrich %>% filter(Type == 'indels'), aes(x=ZT, y=Riskratio)) + geom_boxplot() + 
  xlab('|Z-score| threshold') + ylab('') + ggtitle('indels') +
  gtex_v8_figure_theme() + geom_hline(yintercept=1, color='darkgrey')

sfig_C = ggplot(single_enrich %>% filter(Type == 'SVs'), aes(x=ZT, y=Riskratio)) + geom_boxplot() + 
  xlab('|Z-score| threshold') + ylab('') + ggtitle('SVs') +
  gtex_v8_figure_theme() + geom_hline(yintercept=1, color='darkgrey')

senrich <- plot_grid(sfig_A, sfig_B, sfig_C, labels = c('A','B', 'C'), ncol=3)

ggsave(senrich, file=paste0(data_dir, 'paper_figures/sfig_single.pdf'), width=7.2, height=4,units="in")



 

