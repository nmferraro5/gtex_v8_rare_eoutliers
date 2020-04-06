library(ggplot2)
library(data.table)
library(dplyr)
library(RColorBrewer)
library(cowplot)
library(ggpubr)

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
dcols = c('#7F5A83', '#0D324D', '#BFCDE0')
names(dcols) = c('aseOutliers', 'sOutliers', 'eOutliers')

load(paste0(data_dir, 'gtexV8.outlier.variant.count.per.individual.v8.update.RData'))
count0.9 = filter(sum_melted, variable == 'Outliers w/ high WS')
count0.9$variable = 'Watershed > 0.5'
load(paste0(data_dir, 'gtexV8.outlier.Z3.variant.count.per.individual.v8.update.ws0.9.RData'))

load(paste0(data_dir, 'gtexV8.outlier.Z3.variant.count.per.individual.v8.update.ws.RData'))
sum_melted$variable = sapply(sum_melted$variable, function(x)
  ifelse(x == 'Outliers w/ RV', 'Outlier RV',
         ifelse(x == 'Coloc Outliers w/ RV', 'Coloc Outlier RV',
                ifelse(x == 'Outliers w/ high WS > 0.5', 'Watershed > 0.5',
                       ifelse(x == 'Outliers', 'Outliers', 
                              ifelse(x == 'Outliers w/ high WS > 0.7', 'Watershed > 0.7', 'Watershed > 0.9'))))))


#sum_melted = rbind(sum_melted, count0.9)
sum_melted$variable = factor(sum_melted$variable, levels=c('Outliers', 'Outlier RV', 'Coloc Outlier RV', 'Watershed > 0.5', 'Watershed > 0.7', 'Watershed > 0.9'))
sum_melted$Type = sapply(sum_melted$Type, function(x) ifelse(x == 'ASE', 'aseOutliers',
                                                             ifelse(x == 'TE', 'eOutliers', 'sOutliers')))
sum_melted$Type = factor(sum_melted$Type, levels=c('eOutliers', 'aseOutliers', 'sOutliers'))
sum_melted = filter(sum_melted, variable != 'Coloc Outlier RV', variable != 'Watershed > 0.7')
fig_5A = ggplot(sum_melted, aes(x=variable,y=value+1)) + 
  geom_boxplot(aes(fill=Type),outlier.size = 0.75) + theme_bw() + # removed outlier.alpha=0.5
  xlab('') + ylab('Number per individual') +
  scale_fill_manual(values=dcols) +
  scale_y_log10() +
  gtex_v8_figure_theme() +
  theme(legend.position=c(0.9,0.85),
        legend.title=element_blank(),
        legend.key.height=unit(0.05,'in'),
        axis.text.x=element_text(hjust=1,angle=35))

### panel B ###
ukbb_any_coloc = fread(paste0(data_dir, 'coloc/gwas/gtex.any.outlier.coloc.variant.ukbb.effect.sizes.v8.update.txt'))
my_comparisons <- list(c("Coloc Outlier", "Outlier"), c("Coloc Outlier", "Non-Outlier"),c("Outlier", "Non-Outlier"))
ukbb_any_coloc$Cat = factor(ukbb_any_coloc$Cat, levels=c('Non-Outlier', 'Outlier', 'Coloc Outlier'))
fig_5B = ggplot(ukbb_any_coloc, aes(x=Cat,y=VPercentile)) + 
  geom_violin(fill='#7EB09B',alpha=0.5) + geom_boxplot(width=0.3) + theme_bw() +
  xlab('') + ylab('Variant effect size percentile') +
  stat_compare_means(method="wilcox.test", comparisons=my_comparisons, 
                     method.args=list(alternative="greater"), label="p.signif", tip.length=0) +
  gtex_v8_figure_theme() + theme(axis.text.x=element_text(hjust=1,angle=35))

### panel C ###
all_stats = fread(paste0(data_dir, 'coloc/gwas/coloc_cadd_watershed_lm_scaled.txt'))
fig_5C = ggplot(all_stats, aes(x=Feature, y= abs(Estimate))) + geom_point(size=2) + 
  geom_errorbar(aes(ymin = abs(Estimate) - 1.95*SE, ymax = abs(Estimate) + 1.95*SE), width=0) + 
  theme_bw() + geom_hline(yintercept=0, color='grey') + xlab('') + ylab('') +
  ylab('Effect size') +
  annotate("text", x=1, y=1.75, label='*', cex=4) +
  annotate("text", x=2, y=2.5, label='***', cex=4) +
  gtex_v8_figure_theme() +
  theme(axis.text.x=element_text(hjust=1,angle=35))

### panel D ###
#hit_random_table = fread(paste0(data_dir, 'coloc/gwas/gtexV8.ukbb.watershed.proportion.hits.txt'))
hit_random_table = fread(paste0(data_dir, 'coloc/gwas/gtexV8.ukbb.watershed.proportion.hits.v8.update.txt'))
hit_random_table$PT = factor(hit_random_table$PT, levels=c(0.01, 0.05, 0.25, 0.5, 0.75, 0.9))
fig_5D = ggplot(hit_random_table %>% filter(Cat == 'Random', Model == 'Watershed'), aes(x=PT,y=Prop)) + geom_boxplot() +
  geom_point(data = filter(hit_random_table, Cat == 'Actual', Model == 'Watershed'), aes(x=PT, y=Prop),color = 'red',size=1) +
  theme_bw() + xlab('Posterior threshold') + ylab('Proportion in top 25%') +
  gtex_v8_figure_theme() +
  theme(legend.position=c(0.75,0.85),
        legend.title=element_blank(),
        legend.key.height=unit(0.05,'in'))

### panel E - example 1 ###
cgene = 'ENSG00000137033.11'
gstart = 6215786
gend = 6257983

all_trait_data = fread(paste0(data_dir, 'coloc/gwas/traits/asthma.chr9.ukbb.sum.stats.tsv'))
all_trait_data$Chr = sapply(all_trait_data$variant, function(x) strsplit(x, ':')[[1]][1])
all_trait_data$Pos = sapply(all_trait_data$variant, function(x) strsplit(x, ':')[[1]][2])
all_trait_data$VID = paste(paste0('chr', all_trait_data$Chr), all_trait_data$Pos, sep=':')
mp_data = all_trait_data %>% filter(Chr == 9, low_confidence_variant == 'FALSE') %>% 
  select(Chr, Pos, pval, variant) %>% filter(!is.na(pval), !is.na(Chr))
mp_data$Pos = as.numeric(mp_data$Pos)
mp_data$Chr = as.numeric(mp_data$Chr)
## manhattan plot
mp_data$IsVar = sapply(mp_data$variant, function(x) ifelse(x == '9:6252690:T:G' | x == '9:6255967:G:C', 'VOI', 
                                                           ifelse(x == '9:6210099:T:C', 'leadSNP', 'Other')))
mp_data$IsVar = factor(mp_data$IsVar, levels=c('VOI', 'leadSNP', 'Other'))
fig_5E1 = ggplot(mp_data, aes(x=Pos,y=-log10(pval))) + geom_point(size=0.5) + 
  theme_bw() + guides(color=F) + 
  ggtitle('Asthma') + xlab('Chromosome 9 position') +
  geom_point(data = subset(mp_data, IsVar == 'VOI'),
             aes(x = Pos, y = -log10(pval)), color = '#CC00B4', size = 1) +
  geom_point(data = subset(mp_data, IsVar == 'leadSNP'),
             aes(x = Pos, y = -log10(pval)), color = '#2699FF', size = 1) +
  geom_hline(yintercept=-log10(5e-08), color='grey') +
  gtex_v8_figure_theme()

all_trait_plot_data = filter(all_trait_data, Chr == 9, Pos > gstart - 10000, Pos < gend + 10000)
all_trait_plot_data$VOI = sapply(all_trait_plot_data$VID, function(x) ifelse(x == 'chr9:6252690' | x == 'chr9:6255967', 'Watershed variant',
                                                                             ifelse(x == 'chr9:6210099', 'lead SNP', 'Other')))
all_trait_plot_data$VOI = factor(all_trait_plot_data$VOI, levels=c('Watershed variant','lead SNP','Other'))
cc_ratio = 41633/318894
scale_factor = cc_ratio * (1-cc_ratio)
all_trait_plot_data$beta_scaled = all_trait_plot_data$beta / scale_factor
fig_5E2 = ggplot(all_trait_plot_data %>% filter(low_confidence_variant == 'FALSE'), aes(x=minor_AF,y=abs(beta_scaled))) +
  geom_point(aes(color=VOI),size=0.5) + theme_bw() + guides(size=F) +
  geom_point(data = subset(all_trait_plot_data, VOI != 'Other'),
             aes(color = VOI), size=1) +
  ylab('|Effect size|') +
  scale_color_manual(values=c('#CC00B4', '#2699FF', 'black')) + xlab('UKBB MAF') +
  gtex_v8_figure_theme() + guides(color=F)

fig_5E <- plot_grid(fig_5E1, fig_5E2, labels = NULL, nrow=2)

### panel F - example 2 ###
cgene = 'ENSG00000177700.5'
gstart = 837356
gend = 842545	

all_trait_data = fread(paste0(data_dir, 'coloc/gwas/traits/bmi.chr11.ukbb.sum.stats.tsv'))
all_trait_data$Chr = sapply(all_trait_data$variant, function(x) strsplit(x, ':')[[1]][1])
all_trait_data$Pos = sapply(all_trait_data$variant, function(x) strsplit(x, ':')[[1]][2])
all_trait_data$VID = paste(paste0('chr', all_trait_data$Chr), all_trait_data$Pos, sep=':')
mp_data = all_trait_data %>% filter(Chr == 11, low_confidence_variant == 'FALSE') %>% 
  select(Chr, Pos, pval, variant) %>% filter(!is.na(pval), !is.na(Chr))
mp_data$Pos = as.numeric(mp_data$Pos)
mp_data$Chr = as.numeric(mp_data$Chr)
## manhattan plot
# mp_data$IsVar = sapply(mp_data$variant, function(x) ifelse(x %in% c('11:825341:C:T', '11:840224:T:C', '11:846407:C:T'), 'VOI', 
#                                                            ifelse(x == '11:840363:C:T', 'leadSNP', 'Other')))
mp_data$IsVar = sapply(mp_data$variant, function(x) ifelse(x == '11:840224:T:C', 'VOI', 
                                                           ifelse(x == '11:840363:C:T', 'leadSNP', 'Other')))
mp_data$IsVar = factor(mp_data$IsVar, levels=c('VOI', 'leadSNP', 'Other'))
fig_5F1 = ggplot(mp_data, aes(x=Pos,y=-log10(pval))) + geom_point(size=0.5) + 
  theme_bw() + guides(color=F) + 
  ggtitle('BMI') + xlab('Chromosome 11 position') +
  geom_point(data = subset(mp_data, IsVar == 'VOI'),
             aes(x = Pos, y = -log10(pval)), color = '#CC00B4', size = 1) +
  geom_point(data = subset(mp_data, IsVar == 'leadSNP'),
             aes(x = Pos, y = -log10(pval)), color = '#2699FF', size = 1) +
  geom_hline(yintercept=-log10(5e-08), color='grey') +
  gtex_v8_figure_theme()

all_trait_plot_data = filter(all_trait_data, Chr == 11, Pos > gstart - 15000, Pos < gend + 15000)
# all_trait_plot_data$VOI = sapply(all_trait_plot_data$VID, function(x) ifelse(x %in% c('chr11:825341', 'chr11:840224', 'chr11:846407'), 'Watershed variant',
#                                                                              ifelse(x == 'chr11:840363', 'lead SNP', 'Other')))
all_trait_plot_data$VOI = sapply(all_trait_plot_data$VID, function(x) ifelse(x == 'chr11:840224', 'Watershed variant',
                                                                             ifelse(x == 'chr11:840363', 'lead SNP', 'Other')))

all_trait_plot_data$VOI = factor(all_trait_plot_data$VOI, levels=c('Watershed variant','lead SNP','Other'))
fig_5F2 = ggplot(all_trait_plot_data %>% filter(low_confidence_variant == 'FALSE'), aes(x=minor_AF,y=abs(beta))) +
  geom_point(aes(color=VOI),size=0.4) + theme_bw() + guides(size=F) +
  geom_point(data = subset(all_trait_plot_data, VOI != 'Other'),
             aes(color = VOI),size=1) +
  scale_color_manual(values=c('#CC00B4', '#2699FF', 'black')) + xlab('UKBB MAF') +
  gtex_v8_figure_theme() + ylab('|Effect size|') +
  theme(legend.title=element_blank(),
        legend.position = c(0.75, 0.75),
        legend.key.height=unit(0.05,'in'))

fig_5F <- plot_grid(fig_5F1, fig_5F2, labels = NULL, nrow=2)

## alternate second example
cgene = 'ENSG00000075234.16' # coloc in muscle + fibroblasts, chr22_46287079
gstart = 46267961
gend = 46294008

all_trait_data = fread(paste0(data_dir, 'coloc/gwas/traits/cholesterol.chr22.ukbb.sum.stats.tsv'))
all_trait_data$Chr = sapply(all_trait_data$variant, function(x) strsplit(x, ':')[[1]][1])
all_trait_data$Pos = sapply(all_trait_data$variant, function(x) strsplit(x, ':')[[1]][2])
all_trait_data$VID = paste(paste0('chr', all_trait_data$Chr), all_trait_data$Pos, sep=':')

all_trait_data$Pos = as.numeric(all_trait_data$Pos)
all_trait_plot_data = filter(all_trait_data, Chr == 22, Pos > gstart - 10000, Pos < gend + 10000)
## choose 22:46287046:C:G as surrogate for eQTL signal
all_trait_plot_data$VOI = sapply(all_trait_plot_data$VID, function(x) ifelse(x == 'chr22:46291947', 'Watershed variant','Other'))
all_trait_plot_data$VOI = factor(all_trait_plot_data$VOI, levels=c('Watershed variant', 'eQTL', 'Other'))
cc_ratio = 43957/317184
scale_factor = cc_ratio * (1-cc_ratio)
all_trait_plot_data$beta_scaled = all_trait_plot_data$beta / scale_factor

kind = which(all_trait_data$VID == 'chr22:46291947')
all_trait_data$IsVar = 'Other'
all_trait_data$IsVar[kind] = 'VOI'
fig_5F1 = ggplot(all_trait_data %>% filter(low_confidence_variant == 'FALSE'), aes(x=Pos,y=-log10(pval))) + geom_point(size=0.5) + 
  theme_bw() + guides(color=F) + 
  ggtitle('High Cholesterol') + xlab('Chromosome 22 position') +
  geom_point(data = subset(all_trait_data, IsVar == 'VOI'),
             aes(x = Pos, y = -log10(pval)), color = '#CC00B4', size = 1) +
  geom_point(data = subset(all_trait_data, IsVar == 'leadSNP'),
             aes(x = Pos, y = -log10(pval)), color = '#2699FF', size = 1) +
  geom_hline(yintercept=-log10(5e-08), color='grey') +
  gtex_v8_figure_theme()

fig_5F2 = ggplot(all_trait_plot_data %>% filter(low_confidence_variant == 'FALSE'), aes(x=minor_AF,y=abs(beta_scaled))) +
  geom_point(aes(color=VOI),size=1) + theme_bw() + guides(size=F) +
  geom_point(data = subset(all_trait_plot_data, VOI != 'Other'),
             aes(color = VOI), size=1) +
  ylab('|Effect size|') +
  scale_color_manual(values=c('#CC00B4', 'black')) + xlab('UKBB MAF') +
  gtex_v8_figure_theme() + guides(color=F)

fig_5F <- plot_grid(fig_5F1, fig_5F2, labels = NULL, nrow=2)

## compile
first_row <- plot_grid(fig_5A, fig_5B, fig_5D, labels = c('A','B','C'), ncol=3, axis='tlbr', align='hv')
second_row1 <- plot_grid(fig_5E, labels = c('D'), ncol=1)
second_row2 <- plot_grid(fig_5F, labels = c('E'), ncol=1)
second_row <- plot_grid(second_row1, second_row2, ncol=2)

fig_5 = plot_grid(first_row, NULL, nrow = 2)

ggsave(fig_5, file=paste0(data_dir, 'paper_figures/current_figures/fig5_revisions_v2.svg'), width=7.2, height=6, units="in")
ggsave(second_row, file=paste0(data_dir, 'paper_figures/fig5_revisions_secondrow.png'), width=7.2, height=3, units='in')

#### EXTRAS
### tissue specificity ###
all_tissue_data = fread(paste0(data_dir, 'coloc/gwas/gtex.ukbb.tissue.specific.all.types.outlier.trait.betas.txt'))
all_tissue_data$IsOutlier = factor(all_tissue_data$IsOutlier, levels=c(0,1))
colnames(all_tissue_data)[7] = 'Code'
trait_map = fread(paste0(data_dir, 'coloc/gwas/gtex.gwas.tissue.specific.traits.txt'))
trait_map$Name[4] = 'Hypothyroidism'
all_tissue_data = inner_join(all_tissue_data, trait_map, by='Code')
all_tissue_data$IsOutlier = sapply(all_tissue_data$IsOutlier, function(x)
  ifelse(x == 1, 'Outlier gene', 'Non-outlier gene'))
fig_5F = ggplot(all_tissue_data, aes(x=Name, y=abs(beta))) + geom_boxplot(aes(fill=IsOutlier)) +
  theme_bw() + scale_y_log10() + xlab('') +
  scale_fill_manual(values=c('grey', 'darkcyan')) +
  annotate("text", x=1:13, y=0.05, label='****',cex=4) +
  gtex_v8_figure_theme() +
  theme(axis.text.x=element_text(hjust=1, angle=25),
        legend.title=element_blank(),
        legend.key.height=unit(0.05,'in'),
        legend.position = c(0.95, 0.95))


