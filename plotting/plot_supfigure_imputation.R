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

exp.design = fread(paste0(data_dir, 'gtex_2017-06-05_v8_design_passed_v2.txt'))


load(paste0(data_dir, 'compare_imputation/GTEXv8_pick_parameter_values_and_compare_methods.RData'))
meth.colors = c(brewer.pal(9, 'YlGnBu')[3:9],'grey')
names(meth.colors) = c('EM','KNN','MEAN','PMD','SOFT','STFZ','MEDZ','RANDOM')
keep_methods = c(em.min,knn.min,pmd.min,soft.min[1],'Random','Mean')
method_table = rbind(c('Random','RANDOM'),
                     c('EM_3', 'EM'),
                     c('Mean', 'MEAN'),
                     c('Knn_200', 'KNN'),
                     c('Soft_20_20_15', 'SOFT'),
                     c('PMD_1_5_1', 'PMD'))
rownames(method_table) = method_table[,1]
recon_keep = filter(reconstruction.long, Method %in% keep_methods)
recon_keep$MethodName = sapply(recon_keep$Method, function(x) method_table[x,2])
recon_keep$MethodName = factor(recon_keep$MethodName,levels=c('EM', 'KNN', 'MEAN', 'PMD', 'SOFT', 'RANDOM'))
sfig_A = ggplot(recon_keep,aes(x=MethodName,y=Error,Group=MethodName)) + 
  geom_boxplot(aes(fill=MethodName),outlier.size=0.5) + theme_bw() +
  scale_fill_manual(values=meth.colors) + xlab('Imputation method') +
  guides(fill=F) + gtex_v8_figure_theme() +
  theme(axis.text.x=element_text(hjust=1,angle=40))


recon.em.mean.knn$Cat = sapply(recon.em.mean.knn$Parameter, function(x) ifelse(x == 200, 1, 0))
recon.em.mean.knn$Cat = factor(recon.em.mean.knn$Cat, levels=c(0,1))
sfig_B = ggplot(recon.em.mean.knn %>% filter(Method == 'Knn', !(Parameter %in% c(5,15,25,35,45,55,65))), aes(x = factor(Parameter), y = Error)) +
  geom_boxplot(aes(fill=Cat),outlier.size=0.5) + xlab('Value of k') +
  scale_y_continuous(limits = c(0, 2)) + gtex_v8_figure_theme() +
  scale_fill_manual(values=c('white', '#E5CEDC')) + guides(fill=F) +
  theme(axis.text.x=element_text(hjust=1,angle=90))


knn_enrich = fread(paste0(data_dir, 'enrichments/enrichments_relativeRisk_10kb_genebody.v8ciseQTLs.globalOutliers.globalGenes.removed.subset.knn.none.txt'))
pts = c(1e-05, 1e-07, 1e-09, 1e-11, 1e-13, 1e-15)
knn_enrich$PT = factor(knn_enrich$PT, levels=pts)
sfig_C = ggplot(knn_enrich, aes(x=PT,y=Riskratio,Group=Method)) + 
  geom_point(size=2,aes(color=Method),position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width=0,position=position_dodge(width=0.5)) +
  theme_bw() + ylab('Relative risk') + 
  scale_color_manual(values=c('#7FCDBB', 'black')) +
  ggtitle('') + theme_bw() + xlab('P-value threshold') +
  geom_hline(yintercept=1, linetype="dashed", colour="darkgrey") +
  gtex_v8_figure_theme() +
  facet_wrap(.~Type,scales='free') +
  theme(axis.text.x=element_text(hjust=1,angle=40)) +
  theme(strip.background = element_blank(),
        strip.text.x=element_text(size=8),
        legend.title=element_blank())

first_row <- plot_grid(sfig_A, sfig_B, sfig_C, labels = c('A','B', 'C'), ncol=3, rel_widths=c(1,1,2), axis='tblr', align='h')

ggsave(first_row, file=paste0(data_dir, 'paper_figures/sfig_impute_v2.pdf'), width=7.2, height=3,units="in")





