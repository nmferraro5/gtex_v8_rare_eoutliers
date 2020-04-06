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
load(paste0(data_dir, 'gtexV8.relative.risks.ase.v8.vg.RData'))
ase_risks = risks
load(paste0(data_dir, 'gtexV8.relative.risks.all.types.RData'))
risks = filter(risks, Type != 'ASE')
risks = rbind(risks, ase_risks)
risks$Cat = factor(risks$Cat, levels=c('no_variant','other_noncoding','TE', 'coding','TSS', 'conserved_noncoding','INV','BND','DEL','CNV','DUP', 'splice_region_variant', 'splice_acceptor_variant','frameshift','splice_donor_variant', 'stop'))
plot_cols['no_variant'] = 'darkgrey'

risk_compare = data.frame(FC = numeric(), Type = character(), Cat = character(), MaxType = character())
for (vcat in unique(risks$Cat)) {
  data = filter(risks, Cat == vcat)
  max_type = filter(risks, Risk == max(data$Risk))$Type
  odata = data %>% mutate(FC = Risk / max(Risk)) %>% select(FC,Type,Cat)
  odata$MaxType = max_type
  risk_compare = rbind(risk_compare, odata)
}

risk_compare = risk_compare %>% arrange(by=FC)
risk_compare$Cat = factor(risk_compare$Cat, levels=unique(risk_compare$Cat))
risk_compare$Type = sapply(risk_compare$Type, function(x) ifelse(x == 'ASE', 'aseOutliers',
                                                                 ifelse(x == 'Splicing', 'sOutliers', 'eOutliers')))
risk_compare$Type = factor(risk_compare$Type, levels=c('eOutliers', 'aseOutliers', 'sOutliers'))
sfig_A = ggplot(risk_compare %>% filter(is.finite(FC)), aes(x=Cat, y=FC)) + 
  geom_point(aes(color=Cat,shape=Type),size=2) + theme_bw() + xlab('') +
  scale_shape_manual(values=c('square', 'circle', 'triangle')) +
  ylab('Enrichment / max enrichment') +
  scale_color_manual(values=plot_cols) + 
  guides(color=F) + geom_hline(yintercept=1,color='grey',linetype='dashed') +
  gtex_v8_figure_theme() + theme(axis.text.x=element_text(hjust=1,angle=45),
                                 legend.title=element_blank(),
                                 legend.position=c(0.9,0.1))


### panel B ###
ase_data = fread(paste0(data_dir, 'absolute_outlier_risk_ase_v8_vg_sv7_allGenes_allVariants_update.txt')) %>%
  filter(Feature != 'Feature')
abs_data = fread(paste0(data_dir, 'absolute_outlier_risk_exp_splice_ase_sv7_allGenes_allVariants_update.txt')) %>%
  filter(Feature != 'Feature')
abs_data$ASE_Risk = ase_data$ASE_Risk
abs_risk = melt(abs_data, id.vars='Feature') 
colnames(abs_risk) = c('Feature','variable','Risk')
abs_risk$Risk = as.numeric(abs_risk$Risk)
abs_risk = abs_risk %>% arrange(by=Risk) %>% mutate(Cat = 'Actual')
random_data = melt(fread(paste0(data_dir, 'absolute_outlier_risk_exp_splice_ase_sv7_allGenes_allVariants_random.txt')) %>%
                     filter(Feature != 'Feature'), id.vars='Feature') %>% mutate(Cat = 'Random')
colnames(random_data)[3] = 'Risk'
random_data$Risk = as.numeric(random_data$Risk)
feat_pvals = data.frame(Feature = character(), variable = character(), Pval = numeric())
for (feat in unique(abs_risk$Feature)) {
  adata = filter(abs_risk, Feature == feat)
  rdata = filter(random_data, Feature == feat)
  epval = length(which(filter(rdata, variable == 'MedZ_Risk')$Risk >= filter(adata, variable == 'MedZ_Risk')$Risk)) / 1000 
  spval = length(which(filter(rdata, variable == 'Splice_Risk')$Risk >= filter(adata, variable == 'Splice_Risk')$Risk)) / 1000 
  apval = length(which(filter(rdata, variable == 'ASE_Risk')$Risk >= filter(adata, variable == 'ASE_Risk')$Risk)) / 1000 
  feat_pvals = rbind(feat_pvals, data.frame(Feature = feat, variable = c('MedZ_Risk', 'Splice_Risk', 'ASE_Risk'), Pval= c(epval,spval,apval)))
}
abs_risk = inner_join(abs_risk, feat_pvals, by=c('Feature', 'variable'))
abs_risk$Feature = factor(abs_risk$Feature,levels=unique(abs_risk$Feature))
abs_risk$variable = sapply(abs_risk$variable, function(x)
  ifelse(x == 'MedZ_Risk', 'Total Expression',
         ifelse(x == 'Splice_Risk', 'Splicing', 'ASE')))
random_summary = random_data %>% group_by(Feature,variable) %>%
  mutate(MeanProp = mean(Risk)) %>% mutate(LowProp = mean(Risk) - 1.95*sd(Risk)) %>%
  mutate(HighProp = mean(Risk) + 1.95*sd(Risk)) %>% sample_n(1) %>% ungroup()
dcols = c('#7F5A83', '#0D324D', '#BFCDE0')
names(dcols) = c('aseOutliers', 'sOutliers', 'eOutliers')
abs_risk$variable = sapply(abs_risk$variable, function(x) ifelse(x == 'ASE', 'aseOutliers',
                                                                 ifelse(x == 'Splicing', 'sOutliers', 'eOutliers')))
abs_risk$variable = factor(abs_risk$variable, levels=c('eOutliers', 'aseOutliers', 'sOutliers'))
sfig_B = ggplot(abs_risk %>% filter(Feature != 'INV', Feature != 'no_variant'),aes(x=Feature,y=Risk,Group=variable)) +
  geom_bar(aes(fill=variable),stat='identity',color='black',position='dodge') +
  scale_y_sqrt() + scale_fill_manual(values=dcols) + theme_bw() + 
  xlab('') + ylab('Proportion variants') + 
  gtex_v8_figure_theme() +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        legend.position=c(0.2,0.8),
        legend.title=element_blank(),
        legend.key.height = unit(0.1, "cm"))


## save figure
fig_1 = plot_grid(sfig_A, sfig_B, nrow = 2, labels=c('A', 'B'), align='v')

ggsave(fig_1, file=paste0(data_dir, 'paper_figures/sfig_enrich_revisions_v2.pdf'), width=7.2, height=7.2,units="in")

