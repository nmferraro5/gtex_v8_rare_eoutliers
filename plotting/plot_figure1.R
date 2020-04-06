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

### panel B ###
si_risks = fread(paste0(data_dir, 'gtexV8.all.data.types.relative.risk.10kb.gnomad.txt')) %>% filter(Type == 'SNPs+indels')
sv_risks = fread(paste0(data_dir, 'gtexV8.all.data.types.relative.risk.10kb.txt')) %>% filter(Type == 'HallLabSV')
si_ase_risks = fread(paste0(data_dir, 'gtexV8.all.data.types.snps.indels.relative.risk.10kb.ase.v8.vg.gnomad.txt')) %>% filter(Type == 'SNPs+indels')
sv_ase_risks = fread(paste0(data_dir, 'gtexV8.all.data.types.relative.risk.10kb.ase.v8.vg.txt')) %>% filter(Type == 'HallLabSV')
si_risks = rbind(si_risks %>% filter(Method != 'ASE'), si_ase_risks)
sv_risks = rbind(sv_risks %>% filter(Method != 'ASE'), sv_ase_risks)

si_risks$MafBin = sapply(si_risks$Maf, function(x)
  ifelse(x == 'MAFnovel', 'novel',
         ifelse(x == 'MAFac1', 'single',
                ifelse(x == 'MAFac2', 'double',
                       ifelse(x == 'MAFgnomad0-1only', 'rare', 'low frequency')))))
sv_risks$MafBin = sapply(sv_risks$Maf, function(x)
  ifelse(x == 'MAF0-single', 'single',
         ifelse(x == 'MAFsingle-double', 'double',
                ifelse(x == 'MAF0-1', 'rare',
                       ifelse(x == 'MAF1-5', 'low frequency')))))

risks = rbind(si_risks, sv_risks)
risks$MafBin = factor(risks$MafBin, levels=c('novel', 'single', 'double', 'rare', 'low frequency'))

dcols = c('#7F5A83', '#0D324D', '#BFCDE0')
names(dcols) = c('aseOutliers', 'sOutliers', 'eOutliers')

risks$Type = sapply(risks$Type, function(x) ifelse(x == 'HallLabSV', 'SVs', 'SNVs+indels'))
risks$Type = factor(risks$Type, levels=c('SVs', 'SNVs+indels'))
risks$Method = sapply(risks$Method, function(x) ifelse(x == 'TotalExpression', 'eOutliers', 
                                                       ifelse(x == 'ASE', 'aseOutliers', 'sOutliers')))
risks$Method = factor(risks$Method, levels=c('eOutliers', 'aseOutliers', 'sOutliers'))
fig_1B = ggplot(risks, aes(x=MafBin,y=Riskratio,Group=Method)) +
  geom_point(size=3, aes(shape=Method),position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width=0,position=position_dodge(width=0.5)) +
  theme_bw() + ylab('Relative risk') + xlab('') +
  geom_hline(yintercept=1, linetype="dashed", colour="darkgrey") +
  gtex_v8_figure_theme() + scale_shape_manual(values=c('square', 'circle', 'triangle')) +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        strip.text.x=element_text(size=8)) +
  theme(legend.title=element_blank(), legend.position=c(0.8,0.9)) +
  theme(legend.key.height = unit(0.05, "cm")) +
  facet_wrap(.~Type,scales='free',nrow=2) + theme(strip.background = element_blank())

### panel C ###
load(paste0(data_dir, 'variant.annotation.barplots.RData'))
load(paste0(data_dir, 'ase_v8_vg_variant_annotation_plot.RData'))
colnames(te_over_data)[17] = 'pval_bin'
colnames(te_under_data)[17] = 'pval_bin'
all_bin_data = rbind(plot_ase_bin %>% select(variant_cat,pval_bin,NumCat) %>% mutate(DataType = 'aseOutliers'),
                     splice_data %>% select(variant_cat,pval_bin,NumCat) %>% mutate(DataType = 'sOutliers'),
                     te_over_data %>% select(variant_cat,pval_bin,NumCat) %>% mutate(DataType = 'eOutliers - Over'),
                     te_under_data %>% select(variant_cat,pval_bin,NumCat) %>% mutate(DataType = 'eOutliers - Under'))
all_bin_data$DataType = factor(all_bin_data$DataType, levels=c('eOutliers - Over', 'eOutliers - Under', 'aseOutliers', 'sOutliers'))
fig_1C = ggplot(all_bin_data, aes(x=pval_bin,y=NumCat,Group=variant_cat)) +
  geom_bar(aes(fill=variant_cat),color='black',stat='identity') + 
  theme_bw() + coord_flip() + ylab('') + xlab('') +
  scale_fill_manual(values=splice_cols) + guides(shape=F) +
  ggtitle('') + xlab('P-value bin') + ylab('Proportion with variant') +
  theme(axis.title = element_text(size=8),
        axis.text = element_text(size=8),
        strip.text.x=element_text(size=8)) +
  gtex_v8_figure_theme() + theme(legend.title=element_blank()) +
  theme(legend.key.width = unit(0.2, "cm"),
        legend.key.height = unit(0.05, 'cm')) +
  facet_wrap(.~DataType,nrow=2) + theme(strip.background = element_blank())

### panel D ###
load(paste0(data_dir, 'gtexV8.relative.risks.ase.v8.vg.RData'))
ase_risks = risks
load(paste0(data_dir, 'gtexV8.relative.risks.all.types.RData'))
risks = filter(risks, Type != 'ASE')
risks = rbind(risks, ase_risks)
risks$Cat = factor(risks$Cat, levels=c('no_variant','other_noncoding','TE', 'coding','TSS', 'conserved_noncoding','INV','BND','DEL','CNV','DUP', 'splice_region_variant', 'splice_acceptor_variant','frameshift','splice_donor_variant', 'stop'))
plot_cols['no_variant'] = 'darkgrey'

risks$Type = sapply(risks$Type, function(x) ifelse(x == 'ASE', 'aseOutliers',
                                                   ifelse(x == 'Splicing', 'sOutliers', 'eOutliers')))
risks$Type = factor(risks$Type, levels=c('eOutliers', 'aseOutliers', 'sOutliers'))
fig_1Da = ggplot(risks, aes(x=Cat,y=Risk,Group=Type)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width=0,position=position_dodge(width=0.5)) +
  geom_point(size=2, aes(shape=Type,color=Cat),position=position_dodge(width=0.5)) +
  theme_bw() + xlab('') + ylab('Relative risk') + ylim(c(0,300)) +
  ggtitle('') + guides(color=F,shape=F) +
  scale_color_manual(values=plot_cols) +
  scale_shape_manual(values=c('square', 'circle', 'triangle')) +
  geom_hline(yintercept=1, linetype="dashed", colour="darkgrey") +
  gtex_v8_figure_theme() +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        legend.title=element_blank())

cat_keep = as.data.frame(table(filter(risks, Risk < 40)$Cat)) %>% filter(Freq == 3)
fig_1Db = ggplot(risks %>% filter(Cat %in% cat_keep$Var1), aes(x=Cat,y=Risk,Group=Type)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width=0,position=position_dodge(width=0.5)) +
  geom_point(size=1, aes(shape=Type,color=Cat),position=position_dodge(width=0.5)) +
  theme_bw() + ylab('') + xlab('') + scale_y_continuous(trans='log2') +
  ggtitle('') + guides(color=F,shape=F) +
  scale_color_manual(values=plot_cols) +
  scale_shape_manual(values=c('square', 'circle', 'triangle')) +
  geom_hline(yintercept=1, linetype="dashed", colour="darkgrey") +
  gtex_v8_figure_theme() +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        legend.title=element_blank(),
        panel.border = element_rect(colour = "black", fill=NA))

fig_1D = fig_1Da + annotation_custom(ggplotGrob(fig_1Db), xmin = 0.1, xmax = 12, 
                                     ymin = 30, ymax = 350)

cat_keep = filter(risks, Risk > 20)$Cat
risks$Label = sapply(risks$Type, function(x) ifelse(x == 'eOutliers', 'Expression outliers',
                                                    ifelse(x == 'aseOutliers', 'ASE outliers',
                                                           'Splicing outliers')))
risks$Label = factor(risks$Label, levels=c('Expression outliers', 'ASE outliers', 'Splicing outliers'))
fig_1D_sum = ggplot(risks %>% filter(Cat == 'no_variant' | Cat %in% cat_keep), aes(x=Cat,y=Risk,Group=Label)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width=0,position=position_dodge(width=0.5)) +
  geom_point(size=3, aes(shape=Label,color=Cat),position=position_dodge(width=0.5)) +
  theme_bw() + xlab('') + ylab('Relative risk') + ylim(c(0,300)) +
  ggtitle('') + guides(color=F) +
  scale_color_manual(values=plot_cols) +
  scale_shape_manual(values=c('square', 'circle', 'triangle')) +
  geom_hline(yintercept=1, linetype="dashed", colour="darkgrey") +
  gtex_v8_figure_theme() +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        legend.title=element_blank(),
        legend.position=c(0.2,0.9))

### panel E ###
ase_data = fread(paste0(data_dir, 'absolute_outlier_risk_ase_v8_vg_sv7_allGenes_allVariants_update.txt')) %>%
  filter(Feature != 'Feature')
abs_data = fread(paste0(data_dir, 'absolute_outlier_risk_exp_splice_ase_sv7_allGenes_allVariants_update.txt')) %>%
  filter(Feature != 'Feature')
abs_data$ASE_Risk = ase_data$ASE_Risk
abs_risk = melt(abs_data, id.vars='Feature') 
colnames(abs_risk) = c('Feature','variable','Risk')
abs_risk$Risk = as.numeric(abs_risk$Risk)
abs_risk = abs_risk %>% arrange(by=Risk) 
abs_risk$Feature = factor(abs_risk$Feature,levels=unique(abs_risk$Feature))
abs_risk$variable = sapply(abs_risk$variable, function(x)
  ifelse(x == 'MedZ_Risk', 'Total Expression',
         ifelse(x == 'Splice_Risk', 'Splicing', 'ASE')))
feats_keep = filter(abs_risk, Risk > 0.005)$Feature
abs_risk$variable = sapply(abs_risk$variable, function(x) ifelse(x == 'ASE', 'aseOutliers',
                                                                 ifelse(x == 'Splicing', 'sOutliers', 'eOutliers')))
abs_risk$variable = factor(abs_risk$variable, levels=c('eOutliers', 'aseOutliers', 'sOutliers'))
fig_1E = ggplot(filter(abs_risk,Feature %in% feats_keep),aes(x=Feature,y=Risk,Group=variable)) +
  geom_bar(aes(fill=variable),stat='identity',color='black',position='dodge') +
  scale_fill_manual(values=dcols) + theme_bw() + 
  xlab('') + ylab('Proportion of variants leading to outlier') + 
  gtex_v8_figure_theme() +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        legend.position=c(0.25,0.8),
        legend.title=element_blank(),
        legend.key.height = unit(0.1, "cm"))

first_row <- plot_grid(NULL, fig_1B, labels = c('A','B'), ncol=2, align='h')
second_row <- plot_grid(fig_1D, fig_1E, labels=c('C', 'D'), ncol=2, align='h')
third_row <- plot_grid(fig_1C, labels='E', ncol=1)
fig_1 = plot_grid(first_row, second_row, third_row, nrow = 3, align='v')
ggsave(fig_1, file=paste0(data_dir, 'paper_figures/fig1_v15.eps'), width=7.2, height=11, units="in")





