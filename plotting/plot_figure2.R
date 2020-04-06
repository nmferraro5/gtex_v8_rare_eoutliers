args = commandArgs(trailingOnly=TRUE)

library(ggplot2)
library(ggthemes)
library(data.table)
library(dplyr)
library(RColorBrewer)
library(cowplot)
library(grid)
library(gridExtra)

gtex_v8_figure_theme <- function() {
  return(theme(plot.title = element_text(face="plain",size=8), text = element_text(size=8),
               axis.text=element_text(size=7), panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(),panel.background = element_blank(), 
               axis.line = element_line(colour = "black"), panel.border=element_blank(),
               legend.text = element_text(size=7), legend.title = element_text(size=8)))
}


load('fig2_ab_data_v2.RData')

### panel A ###
dcols = c('#7F5A83', '#0D324D', '#BFCDE0')
names(dcols) = c('ASE', 'Splicing', 'TE')
dist.data$Type = sapply(dist.data$Type, function(x) ifelse(x == 'SNPs', 'SNVs', as.character(x)))
dist.data$Type = factor(dist.data$Type, levels=c('SNVs', 'indels', 'SVs'))
dist.data$DataType = factor(dist.data$DataType, levels=c('ASE', 'Splicing', 'TE'))
fig_2A1 = ggplot(dist.data %>% filter(Riskratio > 1, Type == 'SNVs'), aes(x=Window, y=Riskratio,Group=DataType)) + 
  geom_point(size=1,aes(color=DataType), position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width=0,position=position_dodge(width=0.5)) +
  theme_bw() + ylab('Relative risk') + xlab('') + geom_line(aes(group=DataType),size=0.5) +
  geom_hline(yintercept=1,color='grey') +
  scale_color_manual(values=dcols) + guides(color=F) +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        strip.text.x=element_text(size=8)) +
  gtex_v8_figure_theme() +
  facet_wrap(Type~.,scales='free', ncol=3) + theme(strip.background = element_blank())

fig_2A2 = ggplot(dist.data %>% filter(Riskratio > 1, Type != 'SNVs'), aes(x=Window, y=Riskratio,Group=DataType)) + 
  geom_point(size=1,aes(color=DataType), position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width=0,position=position_dodge(width=0.5)) +
  theme_bw() + ylab('') + geom_line(aes(group=DataType),size=0.4) +
  xlab('Distance upstream from gene (kb)') + ylab('') +
  scale_y_log10() + geom_hline(yintercept=1,color='grey') +
  scale_color_manual(values=dcols) +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        axis.title.x=element_text(hjust=0),
        strip.text.x=element_text(size=8),
        legend.position=c(0.9,0.85),
        legend.key.size=unit(0.05,'in')) +
  gtex_v8_figure_theme() +
  facet_wrap(Type~.,scales='free', ncol=3) + theme(legend.title=element_blank(),
                                                       strip.background = element_blank())

fig_2A <- plot_grid(fig_2A1, fig_2A2, ncol=2, rel_widths=c(1,2))

### panel B ###
vcols = c(brewer.pal(11,'Spectral')[1:7],'grey',brewer.pal(11,'Spectral')[8:11])
names(vcols) = unique(plot.tss.data$promoter_motif)
plot.tss.data$ExpBin = factor(plot.tss.data$ExpBin, levels=c('Under', 'Control', 'Over'))
plot.tss.data$promoter_motif = factor(plot.tss.data$promoter_motif, 
                                      levels=c('other_motif', 'Cmyc', 'CTCF', 'E2F4', 'E2F6', 'ELF1', 'Gabp', 'Nrf1', 'Nrsf', 'PU1', 'SP1', 'Srf', 'Tr4', 'USF1', 'Yy1'))
fig_2B = ggplot(plot.tss.data, aes(x=ExpBin,y=NumBin)) + 
  geom_bar(aes(fill=promoter_motif),stat='identity',color='black') + theme_bw() +
  scale_fill_manual(values=vcols) + xlab('Outlier bin') + 
  ylab('Proportion with RV in promoter') +
  theme(panel.border = element_blank(),
        legend.key.height = unit(0.05, "in")) +
  gtex_v8_figure_theme() + theme(legend.title=element_blank(),
                                 legend.position=c(0.5,0.6))

first_row <- plot_grid(fig_2A, fig_2B, labels = c('A','B'), ncol=2)
second_row <- plot_grid(NULL, labels = c('C'), ncol=1)
third_row <- plot_grid(NULL,NULL, labels = c('D','E'), ncol=2)
fourth_row <- plot_grid(NULL,NULL, labels = c('F'), ncol=1)
fig_2 = plot_grid(first_row,second_row,third_row,fourth_row, nrow = 4, align='v')

ggsave(fig_2, file='fig2.pdf', width=7.2, height=10,units="in")


