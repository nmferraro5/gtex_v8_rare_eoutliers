### analysis for supplemental figure 1 in section 1 for TE

library(ggplot2)
library(ggthemes)
library(data.table)
library(dplyr)
library(RColorBrewer)
library(cowplot)
library(viridis)
library(epitools)

gtex_v8_figure_theme <- function() {
  return(theme(plot.title = element_text(face="plain",size=8), text = element_text(size=8),
               axis.text=element_text(size=7), panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(),panel.background = element_blank(), 
               axis.line = element_line(colour = "black"), 
               panel.border=element_blank(),
               legend.text = element_text(size=7), legend.title = element_text(size=8)))
}

dcols = c('#7F5A83', '#0D324D', '#BFCDE0')
names(dcols) = c('aseOutliers', 'sOutliers', 'eOutliers')
data_dir = '/Users/nicoleferraro/durga_local/data/goats_data/v8_data/'

medz_outliers = fread(paste0(data_dir, 'gtexV8.outliers.v8ciseQTLs.globalOutliers.removed.medz.txt'))
ase_outliers = fread(paste0(data_dir, 'ASE/median.ad.scores.uncorrected.no.global.outliers.v8.vg.outliers.only.tsv'))
sp_outliers = fread(paste0(data_dir, 'splicing/cross_tissue_covariate_method_none_no_global_outliers_emperical_pvalue_gene_level_allpops_outliers.txt'))

sample_annotations = fread(paste0(data_dir, 'GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt'))
sample_annotations = sample_annotations %>% select(SUBJID,SEX,AGE,RACE)
all_outliers = rbind(medz_outliers %>% mutate(value=MedZ,SUBJID=Ind) %>% select(Gene,SUBJID,value) %>% mutate(Cat = 'Expression'),
                     ase_outliers %>% mutate(SUBJID=variable, Gene=V1) %>% select(Gene,SUBJID,value) %>% mutate(Cat = 'ASE'),
                     sp_outliers %>% mutate(Gene=CLUSTER_ID,SUBJID=variable) %>% select(Gene,SUBJID,value) %>% mutate(Cat='Splicing'))
all_outliers = merge(all_outliers, sample_annotations, by='SUBJID', all.x=T)

all_summary = all_outliers %>% group_by(SUBJID,Cat) %>% mutate(N = length(unique(Gene))) %>%
  sample_n(1) %>% ungroup()

## SFig 1A - number by population
all_summary$POP = sapply(all_summary$RACE, function(x)
  ifelse(x == 3, 'White (n=668)',
         ifelse(x == 2, 'Black (n=78)',
                ifelse(x == 1, 'Asian (n=10)', 'Other/Unknown (n=8)'))))
all_summary$POP = factor(all_summary$POP, levels=c('Black (n=78)', 'White (n=668)', 'Asian (n=10)', 'Other/Unknown (n=8)'))
all_summary$Cat = sapply(all_summary$Cat, function(x) ifelse(x == 'ASE', 'aseOutliers',
                                                             ifelse(x == 'Splicing', 'sOutliers', 'eOutliers')))
all_summary$Cat = factor(all_summary$Cat, levels=c('eOutliers', 'aseOutliers', 'sOutliers'))
sfig_1A = ggplot(all_summary, aes(x=POP,y=N,Group=Cat)) + geom_boxplot(aes(fill=Cat)) +
  xlab('') + ylab('# outliers/ind') + scale_fill_manual(values=dcols) +
  gtex_v8_figure_theme() + theme(legend.position=c(0.9,0.8), 
                                 legend.title=element_blank()) +
  theme(legend.key.height = unit(0.05, "cm")) 

## SFig 1B - number by over/under
sfig_1B = ggplot(medz_outliers, aes(MedZ)) + geom_histogram(bins=40) + gtex_v8_figure_theme() +
  xlab('Median Z-score') + ylab('Number of eOutliers') +
  annotate("text", x=7, y=500, label='n=2215',cex=3) +
  annotate("text", x=-7, y=500, label='n=1509',cex=3)

## SFig 1C - effect of eQTL correction on enrichment
load(paste0(data_dir, 'enrichments/enrichments_relativeRisk_10kb_genebody_Z3_FDR0.01_peerTop25.RData'))
peer25_risks = filter(all_logits_top) %>% mutate(Cat = 'Top 25% PEER factors')
load(paste0(data_dir, 'enrichments/enrichments_relativeRisk_10kb_genebody_Z3_FDR0.01_peerTop50.RData'))
peer50_risks = filter(all_logits_top) %>% mutate(Cat = 'Top 50% PEER factors')
load(paste0(data_dir, 'enrichments/enrichments_relativeRisk_10kb_genebody_Z2_FDR0.01_nocorrection.RData'))
peerNone_risks = filter(all_logits_top) %>% mutate(Cat = 'No correction') %>%
  filter(Maf != '0-1', Maf != '1-5')
load(paste0(data_dir, 'enrichments/enrichments_relativeRisk_10kb_genebody_Z2_FDR0.01_knownOnly.RData'))
peerKnown_risks = filter(all_logits_top) %>% mutate(Cat = 'Known covariates')


load(paste0(data_dir, 'enrichments/enrichmentsMEDZ_relativeRisk_10kb_genebody_Z3_FDR0.01allTissues.v8ciseQTLs.globalOutliers.removed.RData'))
si_risks = filter(all_logits_top) %>% mutate(Cat = 'Full PEER + eQTL corrected')
load(paste0(data_dir, 'enrichments/enrichments_relativeRisk_10kb_genebody_Z3_FDR0.01allTissues.nothing.removed.RData'))
si_nocorrect_risks = filter(all_logits_top) %>% mutate(Cat = 'Full PEER')
both_risks = rbind(si_risks, si_nocorrect_risks,peer25_risks,peer50_risks, peerNone_risks, peerKnown_risks)
both_risks$MafBin = sapply(both_risks$Maf, function(x)
  ifelse(x == 'novel', 'novel',
         ifelse(x == 'ac1', 'single',
                ifelse(x == 'ac2', 'double',
                       ifelse(x == 'gnomad0-1only', 'rare', 'low frequency')))))
both_risks$Cat = factor(both_risks$Cat, levels=unique(both_risks$Cat))
both_risks$MafBin = factor(both_risks$MafBin, levels=c('novel', 'single', 'double', 'rare', 'low frequency'))

## Effect of global outlier removal
global_compare = fread(paste0(data_dir, 'gtexV8.eoutliers.full.correct.no.global.correct.with.rv.status.txt'))
orig_outliers = filter(global_compare, abs(MedZ.y) > 3)
orig_controls = filter(global_compare, abs(MedZ.y) < 3, Gene %in% orig_outliers$Gene)
new_outliers = filter(global_compare, abs(MedZ.x) > 3)
new_controls = filter(global_compare, abs(MedZ.x) < 3, Gene %in% new_outliers$Gene)
orig_si = rbind(c(nrow(filter(orig_controls, HasSI == 0)), nrow(filter(orig_controls, HasSI == 1))),
                c(nrow(filter(orig_outliers, HasSI == 0)), nrow(filter(orig_outliers, HasSI == 1))))
orig_sv = rbind(c(nrow(filter(orig_controls, InSV == 1, HasSV == 0)), nrow(filter(orig_controls, InSV == 1, HasSV == 1))),
                c(nrow(filter(orig_outliers, InSV == 1, HasSV == 0)), nrow(filter(orig_outliers, InSV == 1, HasSV == 1))))
new_si = rbind(c(nrow(filter(new_controls, HasSI == 0)), nrow(filter(new_controls, HasSI == 1))),
               c(nrow(filter(new_outliers, HasSI == 0)), nrow(filter(new_outliers, HasSI == 1))))
new_sv = rbind(c(nrow(filter(new_controls, InSV == 1, HasSV == 0)), nrow(filter(new_controls, InSV == 1, HasSV == 1))),
               c(nrow(filter(new_outliers, InSV == 1, HasSV == 0)), nrow(filter(new_outliers, InSV == 1, HasSV == 1))))

orig_si_rr = epitab(orig_si, method='riskratio')$tab[2,]
orig_sv_rr = epitab(orig_sv, method='riskratio')$tab[2,]
new_si_rr = epitab(new_si, method='riskratio')$tab[2,]
new_sv_rr = epitab(new_sv, method='riskratio')$tab[2,]

plot_risks = as.data.frame(rbind(orig_si_rr[5:8], orig_sv_rr[5:8], new_si_rr[5:8], new_sv_rr[5:8]))
plot_risks$Type = c('SNVs+indels', 'SVs', 'SNVs+indels', 'SVs')
plot_risks$Cat = rep(c('Full PEER + eQTL corrected', 'Full PEER (no global outliers) + eQTL corrected'),each=2)
colnames(plot_risks) = c('Riskratio', 'Lower', 'Upper', 'Pval', 'Type', 'Cat')
both_risks$Type = sapply(both_risks$Type, function(x) ifelse(x == 'HallLabSV', 'SVs', 'SNVs+indels'))
plot_risks = rbind(plot_risks, both_risks %>% filter(MafBin == 'rare', Cat != 'Full PEER + eQTL corrected') %>% select(Riskratio, Lower, Upper, Pval, Type, Cat))
plot_risks = plot_risks %>% arrange(by=Riskratio)
plot_risks$Cat = factor(plot_risks$Cat, levels=unique(plot_risks$Cat))
sfig_1C = ggplot(plot_risks %>% filter(Cat != 'Known covariates'), aes(x=Cat, y=Riskratio)) + 
  geom_point(size=2, position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width=0,position=position_dodge(width=0.5)) +
  theme_bw() + ylab('Relative risk') + xlab('') +
  geom_hline(yintercept=1, linetype="dashed", colour="darkgrey") +
  gtex_v8_figure_theme() +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        strip.text.x=element_text(size=8)) +
  theme(legend.key.height = unit(0.05, "cm")) +
  theme(legend.title=element_blank(), legend.position=c(0.8,0.8)) + 
  facet_wrap(.~Type, scales='free') + theme(strip.background = element_blank())


## SFig 1D - enrichments across thresholds
efiles = dir(data_dir, '*gtexV8.all.data.types.snps.indels.thresholds*')
get_enrich_data <- function(ef) {
  edata = fread(paste0(data_dir, ef))
  pt = as.numeric(strsplit(strsplit(ef, 'thresholds.')[[1]][2], '.10kb')[[1]][1])
  edata$Threshold = pt
  return(edata %>% select(-Feature))
}

all_edata = do.call(rbind, lapply(efiles, function(x) get_enrich_data(x)))
all_edata = filter(all_edata, Threshold != 1e-09)
all_edata$Threshold = factor(all_edata$Threshold, levels=c(0.01,0.001,1e-04,1e-05,1e-06,1e-07,1e-08))
all_edata$Maf = sapply(all_edata$Maf, function(x) ifelse(x == 'MAFnovel', 'novel',
                                                         ifelse(x == 'MAFgnomad1-5', 'low frequency', 'rare')))
all_edata$Maf = factor(all_edata$Maf, levels=c('novel', 'rare', 'low frequency'))
all_edata$Method = sapply(all_edata$Method, function(x) ifelse(x == 'ASE', 'aseOutliers',
                                                               ifelse(x == 'Splicing', 'sOutliers', 'eOutliers')))
all_edata$Method = factor(all_edata$Method, levels=c('eOutliers', 'aseOutliers', 'sOutliers'))
sfig_1D = ggplot(all_edata, aes(x=Threshold, y=Riskratio, Group=Method)) + 
  geom_point(size=2, aes(color=Threshold,shape=Method),position=position_dodge(width=0.5)) +
  scale_shape_manual(values=c('square', 'circle', 'triangle')) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width=0,position=position_dodge(width=0.5)) +
  theme_bw() + ylab('Relative risk') + xlab('Outlier threshold') +
  scale_color_viridis(discrete=T,end=0.8) + guides(color=F) +
  geom_hline(yintercept=1, linetype="dashed", colour="darkgrey") +
  gtex_v8_figure_theme() +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        strip.text.x=element_text(size=8)) +
  theme(legend.position=c(0.9,0.7), legend.title=element_blank()) +
  theme(legend.key.height = unit(0.05, "cm")) +
  facet_wrap(.~Maf) + theme(strip.background = element_blank())

first_row <- plot_grid(sfig_1A, sfig_1B, labels = c('A', 'B'), ncol=2, align='hv', axis='tlbr')
second_row <- plot_grid(sfig_1C, sfig_1D, labels = c('C', 'D'), ncol=2, align='hv', axis='tlbr')
sfig <- plot_grid(first_row, second_row, nrow=2, rel_heights = c(1,1.4))

ggsave(sfig, file=paste0(data_dir, 'paper_figures/sfig1_updated_v2.svg'), width=7.2, height=6, units="in")





