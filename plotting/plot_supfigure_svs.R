library(ggplot2)
library(ggthemes)
library(data.table)
library(dplyr)
library(RColorBrewer)
library(cowplot)
library(epitools)

gtex_v8_figure_theme <- function() {
  return(theme(plot.title = element_text(face="plain",size=8), text = element_text(size=8),
               axis.text=element_text(size=6), panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(),panel.background = element_blank(), 
               axis.line = element_line(colour = "black"), 
               panel.border=element_blank(),
               legend.text = element_text(size=7), legend.title = element_text(size=8)))
}


data_dir = '/Users/nicoleferraro/durga_local/data/goats_data/v8_data/'


calc_risk <- function(data, zt) {
  outliers = filter(data, abs(MedZ) > zt) %>% group_by(Gene,Ind) %>% sample_n(1) %>% ungroup()
  controls = filter(data, abs(MedZ) < zt, Gene %in% outliers$Gene) %>% group_by(Gene,Ind) %>% sample_n(1) %>% ungroup()
  counttable = rbind(c(nrow(filter(controls, is.na(AF))), nrow(filter(controls, !is.na(AF)))),
                     c(nrow(filter(outliers, is.na(AF))), nrow(filter(outliers, !is.na(AF)))))
  rr = epitab(counttable,method="riskratio")$tab
  estims = data.frame(Riskratio = rr[2,5], Lower = rr[2,6], Upper = rr[2,7],
                      Pval = rr[2, 8], Threshold = zt)
  return(estims)
}

coding_data = fread(paste0(data_dir, 'gtexV8_coding_svs_with_zscores_outliers_controls.txt'))
coding_enrich = do.call(rbind, lapply(list(3,4,5,6,7,8), function(x) calc_risk(coding_data, x)))
noncoding_data = fread(paste0(data_dir, 'gtexV8_noncoding_svs_with_zscores_outliers_controls.txt'))
noncoding_enrich = do.call(rbind, lapply(list(3,4,5,6,7,8), function(x) calc_risk(noncoding_data, x)))

both_enrich = rbind(coding_enrich %>% mutate(Category = 'Overlapping gene'),
                    noncoding_enrich %>% mutate(Category = 'Non-coding within 10kb'))
both_enrich$Threshold = factor(both_enrich$Threshold, levels=c(3,4,5,6,7,8))
both_enrich = filter(both_enrich, Threshold != 7, Threshold != 8)
sfig_A = ggplot(both_enrich, aes(x=Threshold,y=Riskratio,Group=Category)) +
  geom_point(size=1.5, aes(color=Category),position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width=0,position=position_dodge(width=0.5)) +
  theme_bw() + ylab('Relative risk') + xlab('eOutlier |Median Z-score| threshold') +
  geom_hline(yintercept=1, linetype="dashed", colour="darkgrey") +
  gtex_v8_figure_theme() + 
  theme(legend.title=element_blank(), legend.position=c(0.325,0.95)) +
  theme(legend.key.height = unit(0.05, "cm"),
        legend.key.width = unit(0.04, "cm"))

both_data = rbind(coding_data %>% mutate(Category = 'Coding'),
                  noncoding_data %>% mutate(Category = 'Non-coding'))
both_data = filter(both_data, !is.na(Len))

gene_summary = both_data %>% filter(abs(MedZ) >= 3) %>%
  group_by(Ind,Chr,Start,End) %>% mutate(NG = length(unique(Gene))) %>%
  mutate(NDir = length(unique(sign(MedZ)))) %>% ungroup()

num_genes_affected = as.data.frame(table(gene_summary %>% group_by(Chr,Start,End,NG) %>% 
                                           sample_n(1) %>% ungroup() %>% select(NG)))

sfig_B = ggplot(num_genes_affected, aes(x=Var1,y=Freq)) + geom_bar(stat='identity') +
  gtex_v8_figure_theme() + xlab('Number of affected genes') + ylab('Number of rare SVs') +
  annotate("text", x=1, y=125, label='118', cex=3.5) +
  annotate("text", x=2, y=33, label='26', cex=3.5) +
  annotate("text", x=3, y=11, label='4', cex=3.5) +
  annotate("text", x=4, y=8, label='1', cex=3.5) + gtex_v8_figure_theme()

## alternate B option - number outlier-associated rare SVs per person (over and under)
under_data = both_data %>% group_by(Ind) %>% mutate(NumSV = length(which(MedZ < -3))) %>% sample_n(1) %>% ungroup() %>%
  mutate(Direction = 'Under eOutliers')
over_data = both_data %>% group_by(Ind) %>% mutate(NumSV = length(which(MedZ > 3))) %>% sample_n(1) %>% ungroup() %>%
  mutate(Direction = 'Over eOutliers')
sum_data = rbind(under_data, over_data)
sfig_B = ggplot(sum_data, aes(x=Direction, y=NumSV)) + geom_jitter() + geom_violin() + theme_bw() +
  ylab('# outlier-associated rare SVs per person') + xlab('Direction') +
  annotate("text", x=1, y=6, label='Mean = 0.19', cex=3) +
  annotate("text", x=2, y=7.8, label='Mean = 0.23', cex=3) +
  gtex_v8_figure_theme()

## check annotations
sv_annotations = fread(paste0(data_dir, 'rare_svs/rare_sv_gene_annotations.txt'))
sv_annotations$Ind = paste0('GTEX-', sv_annotations$indiv_id)
colnames(sv_annotations)[1] = 'Gene'
sv_annotations = filter(sv_annotations, variant_cat != 'conserved_noncoding')
sv_data_both = inner_join(both_data %>% filter(abs(MedZ) > 3), sv_annotations %>% select(-indiv_id), by=c('Gene', 'Ind'))

gene_summary_opposite = filter(gene_summary, NG > 1)
gene_data_opposite = sv_data_both %>% select(Gene,Ind,MedZ,Chr,Start,End,chr_numeric,pos,variant_cat) %>%
  filter(Gene %in% gene_summary_opposite$Gene)
gene_data_both = inner_join(gene_data_opposite, gene_data_opposite, by=c('Ind', 'Chr', 'Start', 'End')) %>%
  filter(Gene.x != Gene.y) %>% filter(variant_cat.x == variant_cat.y, pos.x == pos.y)

all_sv_ids = paste(gene_data_both$Chr, gene_data_both$Start, sep=':')
all_sv_ids = paste(all_sv_ids, gene_data_both$End, sep='-')
liftover_svs = cbind(unique(all_sv_ids), fread(paste0(data_dir, 'rare_svs/hglft_genome_3627f_1e9ed0.bed'), header=F))
gene_data_both$SV_ID = all_sv_ids
colnames(liftover_svs) = c('SV_ID', 'LO_ID')
gene_data_both = inner_join(gene_data_both, liftover_svs, by='SV_ID')
gene_data_both$TruePos = sapply(gene_data_both$LO_ID, function(x) as.numeric(strsplit(strsplit(x, ':')[[1]][2], '-')[[1]][1]) + 1)
gene_data_both = filter(gene_data_both, pos.x == TruePos) %>% select(-TruePos)

vcols = c('#FCBBA1', '#EA5948', '#9E0142')
names(vcols) = c('BND', 'DEL', 'DUP')

missing_genes = c('ENSG00000128050.8','ENSG00000128059.8','ENSG00000157426.13', 'ENSG00000174780.15')
missing_data = filter(gene_summary_opposite, Gene %in% missing_genes) %>% mutate(variant_cat = 'DUP') %>%
  mutate(chr_numeric = 4, pos = 56828561)
missing_data_both = inner_join(missing_data, missing_data, by=c('Ind', 'Chr', 'Start', 'End')) %>%
  filter(Gene.x != Gene.y) %>% filter(variant_cat.x == variant_cat.y, pos.x == pos.y) %>%
  mutate(SV_ID = 'chr4:55962394-56804800', LO_ID = 'chr4:56828561-57670966')
missing_data_both = missing_data_both %>% select(colnames(gene_data_both))

gene_data_both = rbind(gene_data_both, missing_data_both)
## Need to remove duplicates
sv_data_all = gene_data_both %>% filter(Gene.x == 'None')
for (svid in unique(gene_data_both$SV_ID)) {
  sv_data = filter(gene_data_both, SV_ID == svid)
  ngenes = length(unique(c(sv_data$Gene.x, sv_data$Gene.y)))
  if (ngenes == 2) {
    sv_data = sv_data[1,]
    x = 10
  } else if (ngenes == 3) {
    sv_data = sv_data[1:3,]
  } else {
    sv_data = sv_data[c(1,2,3,5,6,9),]
  }
  sv_data_all = rbind(sv_data_all, sv_data)
}

sfig_C = ggplot(sv_data_all, aes(x=MedZ.x,y=MedZ.y)) + geom_point(size=1.5,aes(color=variant_cat.x)) + theme_bw() +
  xlab('|Median Z| of Gene 1') + ylab('|Median Z| of Gene 2') +
  gtex_v8_figure_theme() + scale_color_manual(values=vcols) +
  geom_hline(yintercept=3, color='grey', linetype='dashed') +
  geom_hline(yintercept=-3, color='grey', linetype='dashed') +
  geom_vline(xintercept=3, color='grey', linetype='dashed') +
  geom_vline(xintercept=-3, color='grey', linetype='dashed') +
  theme(legend.title=element_blank()) +
  theme(legend.key.height = unit(0.05, "cm"),
        legend.position=c(0.9,0.5))

## save figure
fig_1 = plot_grid(sfig_A, sfig_B, sfig_C, labels = c('A', 'B', 'C'), nrow = 1, axis='tlbr')

ggsave(fig_1, file=paste0(data_dir, 'paper_figures/sfig_sv_effects.svg'), width=7.2, height=3,units="in")

## Check coloc for SVs intersecting multiple genes
sv_coloc = fread('/Users/nicoleferraro/durga_local/data/goats_data/v8_data/rare_svs/rare_svs_affecting_multiple_genes_with_coloc.txt')
sv_coloc$LeadChr = sapply(sv_coloc$lead_snp, function(x) ifelse(is.na(x), x, strsplit(x, '_')[[1]][1]))
sv_coloc$LeadPos = sapply(sv_coloc$lead_snp, function(x) ifelse(is.na(x), x, strsplit(x, '_')[[1]][2]))

sv_coloc = sv_coloc %>% select(-Category,-Len)



