library(ggplot2)
library(ggthemes)
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


data_dir = '/Users/nicoleferraro/durga_local/data/goats_data/v8_data/'


### panel A ###
si_continuous = fread(paste0(data_dir, 'gtexV8.all.data.types.snps.indels.continuous.10kb.gnomad.txt'))
sv_continuous = fread(paste0(data_dir, 'gtexV8.all.data.types.svs.continuous.10kb.gnomad.txt'))
si_continuous$Maf_Bin = sapply(si_continuous$Maf, function(x)
  ifelse(x == 'MAFnovel', 'novel',
         ifelse(x == 'MAFac1', 'single',
                ifelse(x == 'MAFac2', 'double',
                       ifelse(x == 'MAFgnomad0-1only', 'rare', 'low frequency')))))
sv_continuous$Maf_Bin = sapply(sv_continuous$Maf, function(x)
  ifelse(x == 'MAF0-single', 'single',
         ifelse(x == 'MAFsingle-double', 'double',
                ifelse(x == 'MAF0-1', 'rare', 'low frequency'))))

si_continuous$Method = sapply(si_continuous$Method, function(x)
  ifelse(x == 'ASE', 'aseOutliers',
         ifelse(x == 'Splicing', 'sOutliers', 'eOutliers')))

sv_continuous$Method = sapply(sv_continuous$Method, function(x)
  ifelse(x == 'ASE', 'aseOutliers',
         ifelse(x == 'Splicing', 'sOutliers', 'eOutliers')))

si_continuous$Maf_Bin = factor(si_continuous$Maf_Bin, levels=c('novel', 'single', 'double', 'rare', 'low frequency'))
si_continuous$Method = factor(si_continuous$Method, levels=c('eOutliers', 'aseOutliers', 'sOutliers'))
sfig_A = ggplot(si_continuous, aes(x=Maf_Bin,y=Beta,Group=Method)) +
  geom_point(size=3, aes(shape=Method),position=position_dodge(width=0.5)) +
  scale_shape_manual(values=c('square', 'circle', 'triangle')) +
  geom_errorbar(aes(ymin = Beta - 1.95*SE, ymax = Beta + 1.95*SE), width=0,position=position_dodge(width=0.5)) +
  theme_bw() + ylab('Beta') + xlab('') +
  geom_hline(yintercept=0, linetype="dashed", colour="darkgrey") +
  gtex_v8_figure_theme() + guides(shape=F) +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        strip.text.x=element_text(size=8)) +
  theme(legend.title=element_blank(), legend.position=c(0.2,0.9)) +
  theme(legend.key.height = unit(0.05, "cm"))

sv_continuous$Maf_Bin = factor(sv_continuous$Maf_Bin, levels=c('novel', 'single', 'double', 'rare', 'low frequency'))
sfig_B = ggplot(sv_continuous, aes(x=Maf_Bin,y=Beta,Group=Method)) +
  geom_point(size=3, aes(shape=Method),position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = Beta - 1.95*SE, ymax = Beta + 1.95*SE), width=0,position=position_dodge(width=0.5)) +
  theme_bw() + ylab('Beta') + xlab('') +
  geom_hline(yintercept=0, linetype="dashed", colour="darkgrey") +
  gtex_v8_figure_theme() + guides(shape=F) +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        strip.text.x=element_text(size=8)) +
  theme(legend.title=element_blank(), legend.position=c(0.2,0.9)) +
  theme(legend.key.height = unit(0.05, "cm"))


## save figure
fig_1 = plot_grid(sfig_A, sfig_B, nrow = 1, labels=c('A', 'B'), align='v')

ggsave(fig_1, file=paste0(data_dir, 'paper_figures/sfig_enrich_v2.pdf'), width=7.2, height=7.2,units="in")

