## check variants in TOPMed

library(ggplot2)
library(data.table)
library(dplyr)
library(RColorBrewer)
library(cowplot)
library(stringr)


gtex_v8_figure_theme <- function() {
  return(theme(plot.title = element_text(face="plain",size=8), text = element_text(size=8),
               axis.text=element_text(size=7), panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(),panel.background = element_blank(), 
               panel.border = element_blank(),
               axis.line = element_line(colour = "black"), 
               legend.text = element_text(size=7), legend.title = element_text(size=8)))
}

data_dir = '/Users/nicoleferraro/durga_local/data/goats_data/v8_data/forPej/v8_vg/'

ex_cor_data = fread(paste0(data_dir, 'Adipose_SDg_comparison_V7VSV8_Data.tsv'))
sfig_A <- ggplot(ex_cor_data, aes(x=vgs_v7, y=vgs_v8)) + geom_point() +
  geom_abline(color='red') + gtex_v8_figure_theme() +
  xlab(bquote('sqrt('*V^G*') in v7')) +
  ylab(bquote('sqrt('*V^G*') in v8')) 

sp_data = fread(paste0(data_dir, 'GTEX_V7VSV8_Spearman_Correlation_Distribution_Data.tsv'))
sp_data$Cat = 'Tissue'
sfig_B <- ggplot(sp_data, aes(x=Cat,y=spearman_correlation)) +  
  geom_violin(fill='lightgrey') + geom_boxplot(width=0.5) +
  gtex_v8_figure_theme() + xlab('') + ylab('Per tissue spearman correlation of v7 and v8') +
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())

num_data = melt(fread(paste0(data_dir, 'GTEX_V7VSV8_Number_Of_VGs_Available_Data.tsv')), id.vars='tissue')
num_data$Version = sapply(num_data$variable, function(x) ifelse(x == 'number_of_genes_v7', 'v7', 'v8'))
num_data = num_data %>% arrange(by=-value)
num_data$tissue = factor(num_data$tissue, levels=unique(num_data$tissue))
vcols = c('#A39592','#441A33')
tissue_names = fread('/Users/nicoleferraro/durga_local/data/goats_data/v8_data/gtex_tissue_map_jun2019.txt')
tissue_names$FullName = sapply(tissue_names$tissue, function(x) str_replace_all(x, '_', ' '))
num_data$FullName = sapply(num_data$tissue, function(x)
  filter(tissue_names, abbreviation == x)$FullName)
num_data$FullName = factor(num_data$FullName, levels=rev(unique(num_data$FullName)))
sfig_C <- ggplot(num_data, aes(x=FullName, y=value, Group=Version)) +
  geom_bar(stat='identity',aes(fill=Version),position='dodge',alpha=0.8) + gtex_v8_figure_theme() +
  ylab('Number of genes') + xlab('') +
  scale_fill_manual(values=vcols) +
  theme(axis.text.x=element_text(hjust=1,angle=90),
        legend.position=c(0.05,0.9),
        legend.title=element_blank())


first_row <- plot_grid(sfig_A, sfig_B, labels=c('A','B'), ncol=2, align='hv', axis='tblr', rel_widths = c(2,1))
sfig <- plot_grid(first_row, sfig_C, labels=c('', 'C'), nrow=2)

fig_dir = '/Users/nicoleferraro/durga_local/data/goats_data/v8_data/paper_figures/'
ggsave(sfig, file=paste0(fig_dir, 'sfig_vgs_v4.pdf'), width=7.2, height=8, units="in")



