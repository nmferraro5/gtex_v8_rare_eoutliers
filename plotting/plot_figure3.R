library(ggplot2)
library(ggthemes)
library(data.table)
library(RColorBrewer)
library(cowplot)
library(tidyverse)
library(readr)

gtex_v8_figure_theme <- function() {
  return(theme(plot.title = element_text(face="plain",size=8), text = element_text(size=8),
               axis.text=element_text(size=7), panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(),panel.background = element_blank(), 
               axis.line = element_line(colour = "black"), 
               panel.border=element_blank(),
               legend.text = element_text(size=7), legend.title = element_text(size=8)))
}


data_dir = '/Users/nicoleferraro/durga_local/data/goats_data/v8_data/'

### panel A ###
var_dir = '/Users/nicoleferraro/durga_local/active_projects/GOATs_v2/plotting/fig3A/'
class_colors <- c(
  rgb(121, 92, 129, maxColorValue = 255), 
  rgb(23, 50, 75, maxColorValue = 255), 
  rgb(193, 205, 222, maxColorValue = 255)
)

results = fread(paste0(var_dir, 'variant_enrichment_figure_tables/single_tissue_vaiant_enrichment_v8_Vg.tsv'))
outlier_threshold <- c(0.05, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9)
results$sig <- factor(results$sig, levels = outlier_threshold)

results %<>%
  group_by(MAF) %>%
  mutate(med_ratio = median(ratio)) %>%
  ungroup %>%
  mutate(MAF = factor(MAF, levels = unique(MAF)))

results$outlier_type = sapply(results$outlier_type, function(x)
  ifelse(x == 'TotalExpression', 'eOutliers', 
         ifelse(x == 'ASE', 'aseOutliers', 'sOutliers')))

names(class_colors) <- c("aseOutliers", "sOutliers", "eOutliers")
results$outlier_type = factor(results$outlier_type, levels=c('eOutliers', 'aseOutliers', 'sOutliers'))
fig_3A = ggplot(results, aes(sig, ratio, fill = outlier_type)) +
  geom_boxplot(size = .3, outlier.size = .3)  +
  scale_fill_manual(values = class_colors) + 
  geom_hline(yintercept = 1, linetype = "dashed") +
  xlab("Outlier p-value") + ylab("Relative risk") +
  gtex_v8_figure_theme() +
  theme(legend.title = element_blank(),
        legend.position=c(0.4,0.8),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size=unit(0.05,'in'))

### old panel B ###
#load('/Users/nicoleferraro/durga_local/active_projects/GOATs_v2/plotting/fig3B/heatmaps_updated_v8_vg.RData')
#fig_3B = plot_grid(ASE.OA, AS.OA, TE.OA, slegend, ncol = 4, rel_widths = c(0.75,0.75,0.75,0.15))

### new panel B ###
adir = '/Users/nicoleferraro/durga_local/data/goats_data/v8_data/revisions/fromJonah/'
sharing_med <- read_tsv(paste0(adir, "sharing_med.tsv"))
fig_order <- sharing_med %>% 
  filter(type == "ASE") %>% 
  arrange(desc(frac)) %>% 
  .$tissue_name

sharing_med$tissue_name <- factor(sharing_med$tissue_name, levels = fig_order)
sharing_med$type = sapply(sharing_med$type, function(x)
  ifelse(x == 'ASE', 'aseOutliers',
         ifelse(x == 'Splicing', 'sOutliers', 'eOutliers')))
sharing_med$type = factor(sharing_med$type, levels=c('eOutliers', 'aseOutliers', 'sOutliers'))
fig_3B <- 
  ggplot(sharing_med, aes(tissue_name, frac, fill = type)) +
  geom_bar(stat = "identity", position = "dodge") + 
  ylab("Sharing Fraction") + 
  xlab("") +
  scale_fill_manual(values = class_colors) +
  geom_errorbar(aes(min = X1, max = X2), position = position_dodge(.85), width = .4) +
  gtex_v8_figure_theme() +
  theme(axis.text.x = element_text(angle = 70, hjust = 1),
        legend.title=element_blank(),
        legend.position=c(0.9,0.9),
        legend.key.size=unit(0.05,'in')) 

### panel C ###
dat <- read_tsv(paste0(data_dir, "enrichments/single_tissue_NMD_variant_enrichment_v8_Vg.tsv"))

dat <- 
  dat %>%
  group_by(MAF) %>%
  mutate(med_ratio = median(ratio)) %>%
  ungroup %>%
  arrange(med_ratio) %>%
  mutate(MAF = factor(MAF, levels = rev(unique(MAF))), 
         outlier_type = factor(outlier_type, levels = c("ASE", "AS", "TE"))) 


gtex_v8_figure_theme <- function() {
  return(
    theme(
      plot.title = element_text(face="plain",size=8), 
      text = element_text(size=8),
      axis.text=element_text(size=7), 
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.background = element_blank(), 
      axis.line = element_line(colour = "black"), 
      legend.text = element_text(size=7), 
      legend.title = element_text(size=8)
    )
  )
}

dat$outlier_type = sapply(dat$outlier_type, function(x) ifelse(x == 'TE', 'eOutliers',
                                                               ifelse(x == 'ASE', 'aseOutliers', 'sOutliers')))
dat$outlier_type = factor(dat$outlier_type, levels=c('eOutliers', 'aseOutliers', 'sOutliers'))
fig_3C = ggplot(dat, aes(x = MAF, y = ratio, fill = outlier_type)) +
  geom_boxplot(size = .3, outlier.size = .3) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1), 
        #plot.margin = unit(c(0,.1,0,.1), "cm"), 
        axis.title.x = element_blank()) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_fill_manual(values = class_colors) +
  ylab("Relative risk") + guides(fill=F) +
  scale_y_continuous(trans='log2') +
  gtex_v8_figure_theme() 

### panel D ###
meth.colors = brewer.pal(9, 'YlGnBu')[3:9]
names(meth.colors) = c('EM','KNN','MEAN','PMD','SOFT','STFZ','MEDZ')

outliers.top.methods = fread(paste0(data_dir, 'zscore_proportions_medz_cor.txt'))
meth.fills = rep(NA, length(meth.colors))
names(meth.fills) = names(meth.colors)
meth.fills['MEDZ'] = alpha('#C1CDDE', 0.4)
meth.fills['KNN'] = alpha(meth.colors['KNN'], 0)
meth.colors.filtered = meth.colors[c(2,7)]
meth.fills.filtered = meth.fills[c(2,7)]
outliers.top.methods$Method = sapply(outliers.top.methods$Method, function(x) ifelse(x == 'COR', 'Correlation', x))
names(meth.colors.filtered)[1] = 'Correlation'
names(meth.fills.filtered)[1] = 'Correlation'
meth.colors.filtered['MEDZ'] = '#C1CDDE'

fig_3D = ggplot(outliers.top.methods, aes(x = prop, colour = Method, fill = Method)) +
  geom_density(bw = 0.04, size = 0.75) + ylab('Density') +
  xlab('Proportion tissues with |Z| >= 3') +
  scale_color_manual(values = meth.colors.filtered, breaks = names(meth.colors.filtered)) +
  scale_fill_manual(values=meth.fills.filtered) +
  guides(fill = FALSE, colour = guide_legend(override.aes = list(fill = meth.fills.filtered))) +
  theme(panel.border = element_blank()) +
  gtex_v8_figure_theme() +
  theme(legend.position = c(0.5, 0.85),
        legend.direction = 'vertical',
        legend.title=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.key.size=unit(0.05,'in'))

### panel E ###
enh_risks = fread(paste0(data_dir, 'enhancers/KNN.FDR.enhancer.snpindels.500kb.risks.jun1.txt'))
ase_risks = fread(paste0(data_dir, 'enhancers/ASE.enhancer.500kb.risks.txt'))
sp_risks = fread(paste0(data_dir, 'enhancers/Splicing.enhancer.500kb.risks.txt'))
enh_risks = filter(enh_risks, ZT %in% c(3,5,7))
ase_risks$ZT = c(3,3,5,5,7,7)
sp_risks$ZT = c(3,3,5,5,7,7)
all_risks = rbind(enh_risks %>% select(Riskratio,Lower,Upper,Pval,MAF,Type,ZT) %>% mutate(Method = 'corOutliers'), 
                  ase_risks %>% select(-PT) %>% mutate(Method = 'aseOutliers'), 
                  sp_risks %>% select(-PT) %>% mutate(Method = 'sOutliers'))
all_risks$Type = sapply(all_risks$Type, function(x)
  ifelse(x == 'Any', 'Unmatched', x))
all_risks$Zscore = factor(all_risks$ZT, levels=c(3,5,7))
fig_3E = ggplot(all_risks %>% filter(Method == 'corOutliers'), aes(x=Type,y=Riskratio,Group=Zscore)) +
  geom_point(size=2, position=position_dodge(width=0.5),aes(alpha=Zscore),color=meth.colors[2]) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width=0, position=position_dodge(width=0.5)) +
  theme_bw() + ylab('Relative risk') + xlab('Enhancer type') +
  labs(alpha = "Tissue |Z|") +
  ggtitle('') + scale_alpha_manual(values=c(0.5,0.7,0.9)) +
  geom_hline(yintercept=1, linetype="dashed", colour="darkgrey") +
  theme(panel.border = element_blank()) +
  gtex_v8_figure_theme() +
  theme(legend.position = c(0.9, 0.75),
        legend.direction = 'vertical',
        legend.key.size=unit(0.05,'in'))


first_column <- plot_grid(fig_3B, labels = c('A'), nrow=1)
second_column <- plot_grid(fig_3A, fig_3C, fig_3D, fig_3E, labels = c('B', 'C', 'D', 'E'), ncol=4, axis='tblr', align='h')
fig_3 = plot_grid(first_column, second_column, nrow = 2, align='hv', rel_heights = c(1.5,1))

ggsave(fig_3, file=paste0(data_dir, 'paper_figures/fig3_v4_revisions.svg'), width=7.2, height=7.2,units="in")

#ggsave(fig_3, file=paste0(data_dir, 'paper_figures/fig3_v8_bottomrow.pdf'), width=7.2, height=4,units="in")

#ggsave(first_column, file=paste0(data_dir, 'paper_figures/fig3_v8_toprow.png'), width=7.2, height=2,units="in")


