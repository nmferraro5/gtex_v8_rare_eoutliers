library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)

data_dir = '/Users/nicoleferraro/durga_local/data/goats_data/v8_data/'

### UKBB phase 2 stats, filtered to all GTEx rare variants
ukbb_sum_stats = fread(paste0(data_dir, 'coloc/gwas/gtex.coloc.traits.beta.gwas.overlap.txt'))
ukbb_sum_stats$Chr = sapply(ukbb_sum_stats$variant, function(x) strsplit(x, ':')[[1]][1])
ukbb_sum_stats$Pos = sapply(ukbb_sum_stats$variant, function(x) strsplit(x, ':')[[1]][2])
ukbb_sum_stats$VID = paste(paste0('chr',ukbb_sum_stats$Chr), ukbb_sum_stats$Pos, sep=':')

trait_table = fread(paste0(data_dir, 'coloc/gwas/ukbb.trait.table.txt'))
colnames(trait_table)[1] = 'Trait'
ukbb_sum_stats = inner_join(ukbb_sum_stats, trait_table, by='Trait')
ukbb_sum_stats = filter(ukbb_sum_stats, !is.na(beta))
ukbb_sum_stats = ukbb_sum_stats %>% group_by(Trait) %>% mutate(VRank = rank(abs(beta))) %>% 
  mutate(VPercentile = ntile(abs(beta),100)) %>% ungroup()

## outlier variants
te_variants = fread(paste0(data_dir, 'gtexV8.outliers.v8ciseQTLs.globalOutliers.removed.medz.withraresnps.txt'))
ase_variants = fread(paste0(data_dir, 'ASE/gtexV8.outliers.ase.withraresnps.v8.update.txt'))
sp_variants = fread(paste0(data_dir, 'splicing/gtexV8.outliers.splicing.withraresnps.txt'))
all_variants = unique(c(te_variants$VID, ase_variants$VID, sp_variants$VID))
all_outlier_genes = unique(c(te_variants$Gene, ase_variants$Gene, sp_variants$Gene))
traits_keep = fread(paste0(data_dir, 'coloc/gwas/only.ukbb.traits.coloc.enloc.sqtl.map.txt'))$Trait

ukbb_sum_coloc_stats = ukbb_sum_stats %>% filter(Trait %in% traits_keep) %>%
  mutate(RareTE = ifelse(VID %in% te_variants$VID, 1, 0)) %>%
  mutate(RareASE = ifelse(VID %in% ase_variants$VID, 1, 0)) %>%
  mutate(RareSplice = ifelse(VID %in% sp_variants$VID, 1, 0)) %>%
  mutate(RareAny = ifelse(VID %in% all_variants, 1, 0))

rare_vid_map = fread(paste0(data_dir, 'coloc/gwas/ukbb_rare_SNPs_genes_map.txt'))
ukbb_sum_coloc_stats = merge(ukbb_sum_coloc_stats, rare_vid_map, by='VID', all.x=T)

## GTEx coloc results (from Ram)
eqtl_coloc_results = fread(paste0(data_dir, 'coloc/coloc_sig_results_0.5.txt'))
eqtl_coloc_results$Trait = sapply(eqtl_coloc_results$filename, function(x) strsplit(x, '__PM__')[[1]][1])

eqtl_results = fread(paste0(data_dir, 'coloc/enloc_sig_results_0.5.txt'))
eqtl_results$Trait = sapply(eqtl_results$filename, function(x) strsplit(x, '_w_')[[1]][2])
eqtl_results$Trait = sapply(eqtl_results$Trait, function(x) strsplit(x, '_enloc_output.txt')[[1]][1])
sqtl_results = fread(paste0(data_dir, 'splicing/sqtl_enloc_rcp0.5_collapsed_bygene.txt'))
colnames(sqtl_results) = c('Chr', 'Start', 'End', 'rcp', 'Trait', 'Tissue', 'GeneChr', 'GeneStart', 'GeneEnd', 'Gene')
sqtl_results = sqtl_results %>% select(rcp, Trait, Gene)

ukbb_coloc_map = fread(paste0(data_dir, 'coloc/gwas/only.ukbb.traits.coloc.enloc.sqtl.map.txt'))
ukbb_coloc = filter(ukbb_sum_coloc_stats, Trait %in% ukbb_coloc_map$Trait)
ukbb_coloc = inner_join(ukbb_coloc, ukbb_coloc_map, by=c('Trait', 'Description'))

eqtl_results = eqtl_results %>% select(rcp, Trait, gene_id)
colnames(eqtl_results) = c('ColocProb', 'Enloc_name', 'Gene')
eqtl_coloc_results = eqtl_coloc_results %>% select(p4,Trait,gene)
colnames(eqtl_coloc_results) = c('ColocProb', 'Coloc_name', 'Gene')
colnames(sqtl_results) = c('ColocProb', 'Sqtl_name', 'Gene')

ukbb_eqtl_coloc = merge(ukbb_coloc, eqtl_coloc_results, by=c('Coloc_name', 'Gene'),all.x=T) %>%
  mutate(ColocTrait = Coloc_name) %>% select(-Enloc_name,-Sqtl_name,-Coloc_name) %>%
  mutate(QTL = 'Expression')
ukbb_eqtl_enloc = merge(ukbb_coloc, eqtl_results, by=c('Enloc_name', 'Gene'),all.x=T) %>%
  mutate(ColocTrait = Enloc_name) %>% select(-Enloc_name,-Sqtl_name,-Coloc_name) %>%
  mutate(QTL = 'Expression')
ukbb_sqtl = merge(ukbb_coloc, sqtl_results, by=c('Sqtl_name', 'Gene'),all.x=T) %>%
  mutate(ColocTrait = Sqtl_name) %>% select(-Enloc_name,-Sqtl_name,-Coloc_name) %>%
  mutate(QTL = 'Splicing')
ukbb_coloc_all = rbind(ukbb_eqtl_coloc, ukbb_eqtl_enloc, ukbb_sqtl) %>%
  group_by(VID,Gene,ColocTrait) %>% mutate(IsColoc = ifelse(any(!is.na(ColocProb)), 1, 0)) %>%
  sample_n(1) %>% ungroup() 

ukbb_sum = ukbb_coloc_all %>% filter(low_confidence_variant == 'FALSE') %>%
  group_by(VID,Gene,Trait) %>% sample_n(1) %>% ungroup() %>% 
  select(VID,Gene,Trait,VPercentile,RareTE,RareASE,RareSplice,RareAny,IsColoc)

ukbb_sum_melted = melt(ukbb_sum, id.vars=c('VID','VPercentile','Gene','Trait','IsColoc'))
ukbb_sum_melted$value = factor(ukbb_sum_melted$value, levels=c(0,1))
ukbb_sum_melted$IsColoc = factor(ukbb_sum_melted$IsColoc, levels=c(0,1))
ukbb_outlier_coloc = filter(ukbb_sum_melted, IsColoc == 1, value == 1) %>%
  mutate(Cat = 'Coloc Outlier')
ukbb_outlier = filter(ukbb_sum_melted, IsColoc == 0, value == 1) %>%
  mutate(Cat = 'Outlier')
ukbb_not_coloc = filter(ukbb_sum_melted, value == 0) %>%
  mutate(Cat = 'Non-Outlier')
ukbb_all_coloc_cats = rbind(ukbb_outlier_coloc, ukbb_outlier, ukbb_not_coloc)
ukbb_all_coloc_cats$Cat = factor(ukbb_all_coloc_cats$Cat, levels=c('Non-Outlier', 'Outlier', 'Coloc Outlier'))
ukbb_all_coloc_cats$variable = sapply(ukbb_all_coloc_cats$variable, function(x)
  ifelse(x == 'RareTE', 'TE',
         ifelse(x == 'RareASE', 'ASE',
                ifelse(x == 'RareSplice', 'Splicing', 'Any'))))
ukbb_all_coloc_cats$variable = factor(ukbb_all_coloc_cats$variable, levels=c('ASE', 'Splicing', 'TE', 'Any'))

ccols = c('#ACCBE1', '#7EB09B', '#476A6F')

gkeep = unique(c(te_variants$Gene, ase_variants$Gene, sp_variants$Gene))
ukbb_any_coloc = filter(ukbb_all_coloc_cats, variable == 'Any', Gene %in% gkeep) %>%
  group_by(VID,Gene,Trait,Cat) %>% sample_n(1) %>% ungroup()

my_comparisons <- list(c("Coloc Outlier", "Outlier"), c("Coloc Outlier", "Non-Outlier"),c("Outlier", "Non-Outlier"))
ggplot(ukbb_any_coloc, aes(x=Cat,y=VPercentile)) + 
  geom_violin(fill='#7EB09B',alpha=0.5) + geom_boxplot(width=0.3) + theme_bw() +
  xlab('') + ylab('Variant effect size percentile') +
  stat_compare_means(method="wilcox.test", comparisons=my_comparisons, 
                     label="p.signif", tip.length=0) +
  theme(axis.title=element_text(size=18),
        axis.text=element_text(size=16),
        legend.text=element_text(size=16),
        legend.title=element_blank()) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

####################### SCRATCH ###########################

## split into outlier, high WS, high CADD categories
exp_ws = fread(paste0(data_dir, 'coloc/gwas/expression.gam.watershed.ukbb.variants.txt'))
exp_ws$Gene = sapply(exp_ws$sample_names, function(x) strsplit(x, ':')[[1]][2])
sp_ws = fread(paste0(data_dir, 'coloc/gwas/splicing.gam.watershed.ukbb.variants.txt'))
sp_ws$Gene = sapply(sp_ws$sample_names, function(x) strsplit(x, ':')[[1]][2])
ase_ws = fread(paste0(data_dir, 'coloc/gwas/ase.gam.watershed.ukbb.variants.txt'))
ase_ws$Gene = sapply(ase_ws$sample_names, function(x) strsplit(x, ':')[[1]][2])

exp_ws = inner_join(exp_ws, ukbb_sum_stats %>% select(minor_AF,low_confidence_variant,beta,pval,VRank,VPRank,Trait,Description,VID), by='VID')
sp_ws = inner_join(sp_ws, ukbb_sum_stats %>% select(minor_AF,low_confidence_variant,beta,pval,VRank,VPRank,Trait,Description,VID), by='VID')
ase_ws = inner_join(ase_ws, ukbb_sum_stats %>% select(minor_AF,low_confidence_variant,beta,pval,VRank,VPRank,Trait,Description,VID), by='VID')

ws_all = rbind(exp_ws %>% mutate(outlier_pvalue = total_expression_outlier_pvalue, gam_posterior = total_expression_gam_posterior, ws_posterior = total_expression_watershed_posterior) %>% 
                 select(-total_expression_outlier_pvalue, -total_expression_gam_posterior, -total_expression_watershed_posterior) %>% mutate(Type = 'TE'),
               sp_ws %>% mutate(outlier_pvalue = splicing_outlier_pvalue, gam_posterior = splicing_gam_posterior, ws_posterior = splicing_watershed_posterior) %>% 
                 select(-splicing_outlier_pvalue, -splicing_gam_posterior, -splicing_watershed_posterior) %>% mutate(Type = 'Splicing'),
               ase_ws %>% mutate(outlier_pvalue = ase_outlier_pvalue, gam_posterior = ase_gam_posterior, ws_posterior = ase_watershed_posterior) %>% 
                 select(-ase_outlier_pvalue, -ase_gam_posterior, -ase_watershed_posterior) %>% mutate(Type = 'ASE'))

cadd_scores = fread(paste0(data_dir, 'coloc/gwas/gtexV8.ukbb.rare.cadd.scores.txt'))
cadd_scores = cadd_scores %>% group_by(VID) %>% top_n(1,RawScore) %>% sample_n(1) %>% ungroup()

high_ws_variants = unique(filter(ws_all, ws_posterior > 0.9)$VID)
high_cadd_variants = cadd_scores %>% top_n(length(high_ws_variants), RawScore) %>% select(VID)

co_betas = data.frame(Trait = character(), Ranks = numeric(), Category = character())
for (gtrait in unique(ukbb_sum$Trait)) {
  print(gtrait)
  coloc_data = filter(ukbb_sum, Trait == gtrait)
  coloc_outlier_data = filter(ukbb_sum, Trait == gtrait, RareAny == 1)
  coloc_ws_data = filter(ukbb_sum, Trait == gtrait, VID %in% high_ws_variants)
  coloc_cadd_data = filter(ukbb_sum, Trait == gtrait, VID %in% high_cadd_variants$VID)
  other_data = filter(ukbb_other_data, Trait == gtrait, !(VID %in% coloc_data$VID))
  other_outlier_data = filter(ukbb_other_data, RareAny == 1, Trait == gtrait, !(VID %in% coloc_data$VID))
  co_betas = rbind(co_betas, data.frame(Trait = rep(gtrait,nrow(coloc_data)), Ranks = coloc_data$VRank, Category=rep('Coloc',nrow(coloc_data))))
  co_betas = rbind(co_betas, data.frame(Trait = rep(gtrait,nrow(coloc_outlier_data)), Ranks = coloc_outlier_data$VRank, Category=rep('Coloc_Outlier',nrow(coloc_outlier_data))))
  co_betas = rbind(co_betas, data.frame(Trait = rep(gtrait,nrow(coloc_ws_data)), Ranks = coloc_ws_data$VRank, Category=rep('Coloc_Watershed',nrow(coloc_ws_data))))
  co_betas = rbind(co_betas, data.frame(Trait = rep(gtrait,nrow(coloc_cadd_data)), Ranks = coloc_cadd_data$VRank, Category=rep('Coloc_CADD',nrow(coloc_cadd_data))))
  co_betas = rbind(co_betas, data.frame(Trait = rep(gtrait,nrow(other_data)), Ranks = other_data$VRank, Category=rep('Other',nrow(other_data))))
  co_betas = rbind(co_betas, data.frame(Trait = rep(gtrait,nrow(other_outlier_data)), Ranks = other_outlier_data$VRank, Category=rep('Other_Outlier',nrow(other_outlier_data))))
}

co_betas$Category = factor(co_betas$Category, levels=c('Other', 'Other_Outlier', 'Coloc', 'Coloc_Outlier', 'Coloc_Watershed', 'Coloc_CADD'))
ggplot(co_betas %>% filter(!(Category %in% c('Other', 'Other_Outlier', 'Coloc_Outlier'))), aes(x=Category,y=Ranks)) + 
  geom_boxplot(aes(alpha=Threshold),fill='grey') + theme_bw() +
  xlab('') + ylab('Variant effect size ranks in co-localized regions') +
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=16),
        legend.title = element_text(size=16),
        legend.text = element_text(size=16)) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

wilcox.test(filter(co_betas, Category == 'Coloc_Watershed', Threshold == 0.5)$Ranks, 
            filter(co_betas, Category == 'Coloc')$Ranks, alternative='g')

ukbb_sum = melt(ukbb_sum, id.vars=c('VID','VRank'))
ukbb_sum$value = factor(ukbb_sum$value, levels=c(0,1))
ggplot(ukbb_sum, aes(x=variable,y=VRank)) +
  geom_boxplot(width=0.5,aes(fill=value)) +
  theme_bw() + xlab('') + ylab('Variant ranks by effect size') +
  annotate("text",x=1,y=46500,label='*',cex=5) +
  scale_fill_manual(values=c('lightgrey', 'darkcyan')) +
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=16),
        legend.text=element_text(size=16),
        legend.title=element_blank()) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

for (vcat in unique(ukbb_sum$variable)) {
  print(vcat)
  yesdata = filter(ukbb_sum, variable == vcat, value == 1)$VRank
  nodata = filter(ukbb_sum, variable == vcat, value == 0)$VRank
  print(wilcox.test(yesdata, nodata, alternative='g')$p.value)
}

## outlier variants in coloc regions
te_genes = unique(te_variants$Gene)
ase_genes = unique(ase_variants$Gene)
sp_genes = unique(sp_variants$Gene)
outlier_genes = unique(c(te_genes, ase_genes, sp_genes))
other_genes = fread(paste0(data_dir, 'v8_genes.txt'),header=F) %>% filter(!(V1 %in% outlier_genes))
coloc_genes = unique(c(ukbb_eqtl_coloc$Gene, ukbb_eqtl_enloc$Gene, ukbb_sqtl$Gene))
no_no = length(which(!(other_genes$V1 %in% coloc_genes)))
no_yes = length(which(other_genes$V1 %in% coloc_genes))
yes_no = length(which(!(outlier_genes %in% coloc_genes)))
yes_yes = length(which(outlier_genes %in% coloc_genes))
counttable = rbind(c(no_no, no_yes), c(yes_no, yes_yes))
all_rr = epitab(counttable, method = 'riskratio')$tab


other_genes = fread(paste0(data_dir, 'v8_genes.txt'),header=F) %>% filter(!(V1 %in% ase_genes))
coloc_genes = unique(c(ukbb_eqtl_enloc$Gene))
no_no = length(which(!(other_genes$V1 %in% coloc_genes)))
no_yes = length(which(other_genes$V1 %in% coloc_genes))
yes_no = length(which(!(ase_genes %in% coloc_genes)))
yes_yes = length(which(ase_genes %in% coloc_genes))
counttable = rbind(c(no_no, no_yes), c(yes_no, yes_yes))
ase_rr = epitab(counttable, method = 'riskratio')$tab

combined_rr = rbind(as.data.frame(t(as.data.frame(all_rr[2,]))) %>% mutate(Type = 'All'), as.data.frame(t(as.data.frame(te_rr[2,]))) %>% mutate(Type = 'TE'),as.data.frame(t(as.data.frame(ase_rr[2,]))) %>% mutate(Type = 'ASE'), as.data.frame(t(as.data.frame(sp_rr[2,]))) %>% mutate(Type = 'Splicing'))
ggplot(combined_rr, aes(x=Type,y=riskratio)) + geom_point(size=3) + 
  theme_bw() + geom_errorbar(aes(ymin = lower, ymax = upper), width=0) +
  geom_hline(yintercept=1,color='grey',linetype='dashed') + xlab('') +
  ylab('Relative risk of outlier gene co-localizing') +
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=16)) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
ukbb_coloc_sum = melt(ukbb_coloc_all %>% select(VID,VRank,VPRank,RareTE,RareASE,RareSp,RareAny), id.vars=c('VID','VRank','VPRank'))
ukbb_coloc_sum$value = factor(ukbb_coloc_sum$value, levels=c(0,1))
ggplot(ukbb_coloc_sum, aes(x=variable,y=VRank)) + 
  geom_boxplot(width=0.5,aes(color=value)) +
  theme_bw() + xlab('') + ylab('Variant rank by effect size') +
  annotate("text", x=1, y=46200, label='****', cex=5) +
  annotate("text", x=3, y=46200, label='****', cex=5) +
  annotate("text", x=4, y=46200, label='*', cex=5) +
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=16),
        legend.text=element_text(size=16),
        legend.title=element_blank()) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

for (vcat in unique(ukbb_coloc_sum$variable)) {
  print(vcat)
  yesdata = filter(ukbb_coloc_sum, variable == vcat, value == 1)$VPRank
  nodata = filter(ukbb_coloc_sum, variable == vcat, value == 0)$VPRank
  print(wilcox.test(yesdata, nodata, alternative='g')$p.value)
}

coloc_variants = unique(ukbb_coloc_all$VID)
noncoloc_variants = unique(filter(ws_all, !(VID %in% coloc_variants))$VID)
ws_cadd_risks = data.frame(Riskratio = numeric(), Lower = numeric(), Upper = numeric(), Pval = numeric(), WT = numeric(), Model = character())
wts = c(0.5, 0.6, 0.7, 0.8, 0.9)
for (wt in wts) {
  print(wt)
  high_ws_genes = unique(filter(ws_all, ws_posterior > wt)$Gene)
  high_ws_variants = unique(filter(ws_all, ws_posterior > wt)$VID)
  low_ws_variants = unique(filter(ws_all, !(VID %in% high_ws_variants), !(Gene %in% high_ws_genes))$VID)
  high_cadd_variants = cadd_scores %>% filter(VID %in% ws_all$VID) %>%
    top_n(length(high_ws_variants),RawScore) %>% select(VID)
  high_cadd_genes = filter(ws_all, VID %in% high_cadd_variants)$Gene
  low_cadd_variants = unique(filter(ws_all, !(VID %in% high_cadd_variants$VID), !(Gene %in% high_cadd_genes))$VID)
  no_no_ws = length(which(low_ws_variants %in% noncoloc_variants))
  no_yes_ws = length(which(low_ws_variants %in% coloc_variants))
  yes_no_ws = length(which(high_ws_variants %in% noncoloc_variants))
  yes_yes_ws = length(which(high_ws_variants %in% coloc_variants))
  ws_counttable = rbind(c(no_no_ws, no_yes_ws), c(yes_no_ws, yes_yes_ws))
  ws_rr = epitab(ws_counttable, method='riskratio')$tab
  ws_cadd_risks = rbind(ws_cadd_risks, data.frame(Riskratio = ws_rr[2,5], Lower = ws_rr[2,6], Upper = ws_rr[2,7], Pval = ws_rr[2,8], WT = wt, Model = 'Watershed'))
  
  no_no_cadd = length(which(low_cadd_variants %in% noncoloc_variants))
  no_yes_cadd = length(which(low_cadd_variants %in% coloc_variants))
  yes_no_cadd = length(which(high_cadd_variants$VID %in% noncoloc_variants))
  yes_yes_cadd = length(which(high_cadd_variants$VID %in% coloc_variants))
  cadd_counttable = rbind(c(no_no_cadd, no_yes_cadd), c(yes_no_cadd, yes_yes_cadd))
  cadd_rr = epitab(cadd_counttable, method='riskratio')$tab
  ws_cadd_risks = rbind(ws_cadd_risks, data.frame(Riskratio = cadd_rr[2,5], Lower = cadd_rr[2,6], Upper = cadd_rr[2,7], Pval = cadd_rr[2,8], WT = wt, Model = 'CADD'))
}

ws_cadd_risks$WT = factor(ws_cadd_risks$WT, levels=c(0.5,0.6,0.7,0.8,0.9))
ggplot(ws_cadd_risks, aes(x=WT,y=Riskratio,Group=Model)) + geom_point(size=3,aes(color=Model),position=position_dodge(width=0.5)) + 
  theme_bw() + geom_errorbar(aes(ymin = Lower, ymax = Upper), width=0,position=position_dodge(width=0.5)) +
  geom_hline(yintercept=1,color='grey',linetype='dashed') + xlab('Posterior threshold') +
  ylab('Relative risk of variant in co-localizing region') +
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=16),
        legend.title=element_blank(),
        legend.text=element_text(size=16)) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


watershed_per_variant = ws_all %>% group_by(VID) %>% top_n(1,ws_posterior) %>% 
  sample_n(1) %>% ungroup() %>% select(VID,ws_posterior)
watershed_per_variant = as.data.frame(watershed_per_variant)
rownames(watershed_per_variant) = watershed_per_variant$VID
gam_per_variant = ws_all %>% group_by(VID) %>% top_n(1,gam_posterior) %>% 
  sample_n(1) %>% ungroup() %>% select(VID,gam_posterior)
gam_per_variant = as.data.frame(gam_per_variant)
rownames(gam_per_variant) = gam_per_variant$VID
cadd_scores = as.data.frame(cadd_scores) %>% select(RawScore,VID) %>% distinct() %>%
  group_by(VID) %>% top_n(1,RawScore) %>% sample_n(1) %>% ungroup()
rownames(cadd_scores) = cadd_scores$VID
ukbb_coloc_ws = filter(ukbb_coloc_all, VID %in% ws_all$VID)

ukbb_sum_stats$WS = sapply(ukbb_sum_stats$VID, function(x) watershed_per_variant[x,2])
ukbb_sum_stats$GAM = sapply(ukbb_sum_stats$VID, function(x) gam_per_variant[x,2])
ukbb_sum_stats$CADD = sapply(ukbb_sum_stats$VID, function(x) as.numeric(cadd_scores[x,1]))

## effect size distributions overall

outlier_stats = filter(ukbb_sum_stats, low_confidence_variant == 'FALSE', RareAny == 1) %>% group_by(VID) %>% top_n(1,VRank) %>% sample_n(1) %>% ungroup()
high_ws = filter(ukbb_sum_stats, low_confidence_variant == 'FALSE', WS > 0.9) %>% group_by(VID) %>% top_n(1,VRank) %>% sample_n(1) %>% ungroup()
high_gam = filter(ukbb_sum_stats, low_confidence_variant == 'FALSE', GAM > 0.135) %>% group_by(VID) %>% top_n(1,VRank) %>% sample_n(1) %>% ungroup()
high_cadd = filter(ukbb_sum_stats, low_confidence_variant == 'FALSE', CADD > 5.3) %>% group_by(VID) %>% top_n(1,VRank) %>% sample_n(1) %>% ungroup()
remaining = filter(ukbb_sum_stats, low_confidence_variant == 'FALSE', !(VID %in% high_ws$VID), 
                   !(VID %in% high_gam$VID), !(VID %in% high_cadd$VID), 
                   !(VID %in% outlier_stats$VID)) %>% group_by(VID) %>% 
  top_n(1,VRank) %>% sample_n(1) %>% ungroup()
vars_plot = rbind(outlier_stats %>% select(VRank) %>% mutate(Cat = 'Outlier'),
                  high_ws %>% select(VRank) %>% mutate(Cat = 'Watershed'),
                  high_gam %>% select(VRank) %>% mutate(Cat = 'GAM'),
                  high_cadd %>% select(VRank) %>% mutate(Cat = 'CADD'),
                  remaining %>% select(VRank) %>% mutate(Cat = 'Other'))
vars_plot$Cat = factor(vars_plot$Cat, levels=c('Other', 'CADD', 'Outlier', 'GAM', 'Watershed'))
ggplot(vars_plot, aes(x=Cat,y=VRank)) + geom_boxplot() + theme_bw() +
  xlab('') + ylab('Best effect size rank per variant') +
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=16),
        legend.title=element_blank(),
        legend.text=element_text(size=16)) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

## now just in co-localized regions
ukbb_coloc_sum = inner_join(ukbb_sum_stats %>% select(VID,Trait,RareAny,WS,GAM,CADD,VRank,VPRank), ukbb_coloc_all %>% select(-RareAny,-VRank,-VPRank), by=c('VID', 'Trait')) %>%
  filter(low_confidence_variant == 'FALSE') %>% group_by(VID,Trait) %>% sample_n(1) %>% ungroup()

outlier_stats = filter(ukbb_coloc_sum, RareAny == 1)
high_ws = filter(ukbb_coloc_sum, WS > 0.9)
high_gam = ukbb_coloc_sum %>% top_n(nrow(high_ws),GAM)
high_cadd = ukbb_coloc_sum %>% top_n(nrow(high_ws), CADD)
remaining = filter(ukbb_coloc_sum, !(VID %in% high_ws$VID), !(VID %in% high_gam$VID), !(VID %in% high_cadd$VID), !(VID %in% outlier_stats$VID))
vars_plot = rbind(outlier_stats %>% select(VRank) %>% mutate(Cat = 'Outlier'),
                  high_ws %>% select(VRank) %>% mutate(Cat = 'Watershed'),
                  high_gam %>% select(VRank) %>% mutate(Cat = 'GAM'),
                  high_cadd %>% select(VRank) %>% mutate(Cat = 'CADD'),
                  remaining %>% select(VRank) %>% mutate(Cat = 'Other'))
vars_plot$Cat = factor(vars_plot$Cat, levels=c('Other', 'CADD', 'Outlier', 'GAM', 'Watershed'))
ggplot(vars_plot, aes(x=Cat,y=VRank)) + geom_boxplot() + theme_bw() +
  xlab('') + ylab('Variant effect size rank in co-localized regions') +
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=16),
        legend.title=element_blank(),
        legend.text=element_text(size=16)) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

get_vranks <- function(gtrait) {
  print(gtrait)
  coloc_trait_data = filter(ukbb_coloc_ws, Trait == gtrait) %>% group_by(VID) %>% sample_n(1) %>% ungroup()
  coloc_trait_data$WS = sapply(coloc_trait_data$VID, function(x) watershed_per_variant[x,2])
  coloc_trait_data$GAM = sapply(coloc_trait_data$VID, function(x) gam_per_variant[x,2])
  coloc_trait_data$CADD = sapply(coloc_trait_data$VID, function(x) as.numeric(cadd_scores[x,1]))
  best_ws_rank = coloc_trait_data %>% top_n(1,WS) %>% sample_n(1) 
  best_cadd_rank = coloc_trait_data %>% top_n(1,CADD) %>% sample_n(1)
  best_gam_rank = coloc_trait_data %>% top_n(1,GAM) %>% sample_n(1)
  best_outlier_rank = coloc_trait_data %>% filter(VID %in% all_variants) %>% top_n(1,VRank) %>% sample_n(1)
  if (nrow(best_outlier_rank) > 0) {
    return(as.data.frame(rbind(cbind(best_ws_rank$VPRank, best_ws_rank$WS, best_ws_rank$minor_AF, 'Watershed', gtrait),
                               cbind(best_cadd_rank$VPRank, best_cadd_rank$CADD, best_cadd_rank$minor_AF, 'CADD', gtrait),
                               cbind(best_gam_rank$VPRank, best_gam_rank$GAM, best_gam_rank$minor_AF, 'GAM', gtrait),
                               cbind(best_outlier_rank$VPRank, NA, best_outlier_rank$minor_AF, 'Outlier', gtrait))))
  }
}

all_coloc_vranks = do.call(rbind, lapply(unique(ukbb_coloc_ws$Trait), get_vranks))
colnames(all_coloc_vranks) = c('VRank', 'Value', 'Minor_AF', 'Model','Trait')
all_coloc_vranks$VRank = as.numeric(as.character(all_coloc_vranks$VRank))
all_coloc_vranks$Minor_AF = as.numeric(as.character(all_coloc_vranks$Minor_AF))
all_coloc_vranks$Value = as.numeric(as.character(all_coloc_vranks$Value))

ggplot(all_coloc_vranks, aes(x=Model, y=VRank)) + geom_boxplot() + theme_bw() +
  xlab('') + ylab('Effect size rank of best variant per model per trait') +
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=16)) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


get_vranks <- function(gtrait) {
  print(gtrait)
  coloc_trait_data = filter(ukbb_coloc_all, low_confidence_variant == 'FALSE', Trait == gtrait) %>% group_by(VID) %>% sample_n(1) %>% ungroup()
  #coloc_trait_data$WS = sapply(coloc_trait_data$VID, function(x) watershed_per_variant[x,2])
  #coloc_trait_data$GAM = sapply(coloc_trait_data$VID, function(x) gam_per_variant[x,2])
  #coloc_trait_data$CADD = sapply(coloc_trait_data$VID, function(x) as.numeric(cadd_scores[x,1]))
  other_trait_data = filter(ukbb_sum_stats, low_confidence_variant == 'FALSE', Trait == gtrait, !(VID %in% coloc_trait_data$VID)) %>% group_by(VID) %>% sample_n(1) %>% ungroup()
  max_other_rank = max(other_trait_data$VRank)
  max_other_prank = max(other_trait_data$VPRank)
  max_coloc_rank = max(coloc_trait_data$VRank)
  max_coloc_prank = max(coloc_trait_data$VPRank)
  return(as.data.frame(rbind(c(max_other_rank, 'EffectSize', 'Other', gtrait),
                             c(max_other_prank, 'Pvalue', 'Other', gtrait),
                             c(max_coloc_rank, 'EffectSize', 'Coloc', gtrait),
                             c(max_coloc_prank, 'Pvalue', 'Coloc', gtrait))))
}

all_coloc_vranks = do.call(rbind, lapply(unique(ukbb_coloc_all$Trait), get_vranks))
colnames(all_coloc_vranks) = c('Rank', 'Type', 'Coloc', 'Trait')
all_coloc_vranks$Rank = as.numeric(as.character(all_coloc_vranks$Rank))
ggplot(all_coloc_vranks, aes(x=Type,y=Rank)) + geom_boxplot(aes(fill=Coloc)) +
  theme_bw()


########################## OLD ########################## 
best_near_gene_rank = filter(ukbb_sum_stats, NearGene == 1, Trait %in% ukbb_coloc_all$Trait) %>% 
  filter(low_confidence_variant == 'FALSE') %>%
  group_by(VID) %>% top_n(1,VRank) %>% sample_n(1) %>% 
  ungroup() %>% mutate(Cat = 'Gene')
best_near_outlier_gene_rank = filter(ukbb_sum_stats, NearOutlierGene == 1, Trait %in% ukbb_coloc_all$Trait) %>%
  filter(low_confidence_variant == 'FALSE') %>%
  group_by(VID) %>% top_n(1,VRank) %>% sample_n(1) %>%
  ungroup() %>% mutate(Cat = 'Outlier Gene')
best_near_coloc_gene_rank = ukbb_coloc_all %>% filter(low_confidence_variant == 'FALSE') %>% 
  group_by(VID) %>% top_n(1,VRank) %>% sample_n(1) %>% ungroup() %>%
  mutate(Cat = 'Coloc Gene')

best_ranks = rbind(best_near_gene_rank %>% select(Cat,VPRank,VRank), 
                   best_near_outlier_gene_rank %>% select(Cat,VPRank,VRank), 
                   best_near_coloc_gene_rank%>% select(Cat,VPRank,VRank)) %>%
  mutate(Model = 'All')

high_ws_variants = filter(ws_all, ws_posterior > 0.5)$VID

best_near_gene_rank = filter(ukbb_sum_stats, NearGene == 1, Trait %in% ukbb_coloc_all$Trait, VID %in% high_ws_variants) %>% 
  filter(low_confidence_variant == 'FALSE') %>%
  group_by(VID) %>% top_n(1,VRank) %>% sample_n(1) %>% 
  ungroup() %>% mutate(Cat = 'Gene')
best_near_outlier_gene_rank = filter(ukbb_sum_stats, NearOutlierGene == 1, Trait %in% ukbb_coloc_all$Trait, VID %in% high_ws_variants) %>%
  filter(low_confidence_variant == 'FALSE') %>%
  group_by(VID) %>% top_n(1,VRank) %>% sample_n(1) %>%
  ungroup() %>% mutate(Cat = 'Outlier Gene')
best_near_coloc_gene_rank = ukbb_coloc_all %>% filter(low_confidence_variant == 'FALSE', VID %in% high_ws_variants) %>% 
  group_by(VID) %>% top_n(1,VRank) %>% sample_n(1) %>% ungroup() %>%
  mutate(Cat = 'Coloc Gene')
best_ws_ranks = rbind(best_near_gene_rank %>% select(Cat,VPRank,VRank), 
                   best_near_outlier_gene_rank %>% select(Cat,VPRank,VRank), 
                   best_near_coloc_gene_rank%>% select(Cat,VPRank,VRank)) %>%
  mutate(Model = 'Watershed')

high_gam_variants = filter(ws_all, gam_posterior > 0.03)$VID

best_near_gene_rank = filter(ukbb_sum_stats, NearGene == 1, Trait %in% ukbb_coloc_all$Trait, VID %in% high_gam_variants) %>% 
  filter(low_confidence_variant == 'FALSE') %>%
  group_by(VID) %>% top_n(1,VRank) %>% sample_n(1) %>% 
  ungroup() %>% mutate(Cat = 'Gene')
best_near_outlier_gene_rank = filter(ukbb_sum_stats, NearOutlierGene == 1, Trait %in% ukbb_coloc_all$Trait, VID %in% high_gam_variants) %>%
  filter(low_confidence_variant == 'FALSE') %>%
  group_by(VID) %>% top_n(1,VRank) %>% sample_n(1) %>%
  ungroup() %>% mutate(Cat = 'Outlier Gene')
best_near_coloc_gene_rank = ukbb_coloc_all %>% filter(low_confidence_variant == 'FALSE', VID %in% high_gam_variants) %>% 
  group_by(VID) %>% top_n(1,VRank) %>% sample_n(1) %>% ungroup() %>%
  mutate(Cat = 'Coloc Gene')
best_gam_ranks = rbind(best_near_gene_rank %>% select(Cat,VPRank,VRank), 
                      best_near_outlier_gene_rank %>% select(Cat,VPRank,VRank), 
                      best_near_coloc_gene_rank%>% select(Cat,VPRank,VRank)) %>%
  mutate(Model = 'GAM')

all_ranks = rbind(best_ws_ranks,best_ranks,best_gam_ranks)
all_ranks$Cat = factor(all_ranks$Cat, levels=c('Gene', 'Outlier Gene', 'Coloc Gene'))
ggplot(all_ranks, aes(x=Cat,y=VRank)) + 
  geom_boxplot(aes(fill=Model)) + theme_bw() + xlab('') +
  ylab('Best effect size rank per variant') +
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=16),
        legend.text=element_text(size=16),
        legend.title=element_blank()) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

high_ws_variants = filter(ws_all, ws_posterior > 0.5)$VID
high_gam_variants = filter(ws_all, gam_posterior > 0.14)$VID # 0.14 to match WS > 0.9, 0.0535 to match > 0.7
high_cadd_variants = filter(cadd_scores, RawScore > 5.975)$VID # 5.975 to match WS > 0.9, >5 to match 0.7

ukbb_coloc_all = rbind(ukbb_eqtl_coloc, ukbb_eqtl_enloc) %>%
  mutate(HighWS = ifelse(VID %in% high_exp_variants, 1, 0)) %>%
  mutate(HighGAM = ifelse(VID %in% high_gam_variants, 1, 0)) %>%
  mutate(HighCADD = ifelse(VID %in% high_cadd_variants, 1, 0)) 

ukbb_coloc_ws = filter(ukbb_coloc_all, HighWS == 1) %>% group_by(VID,Trait) %>% sample_n(1) %>% ungroup()
ukbb_coloc_gam = filter(ukbb_coloc_all, HighGAM == 1) %>% group_by(VID,Trait) %>% sample_n(1) %>% ungroup()
ukbb_coloc_cadd = filter(ukbb_coloc_all, HighCADD == 1) %>% group_by(VID,Trait) %>% sample_n(1) %>% ungroup()

ukbb_coloc_types = rbind(ukbb_coloc_ws %>% select(VPRank) %>% mutate(Cat = 'Watershed'),
                         ukbb_coloc_gam %>% select(VPRank) %>% mutate(Cat = 'GAM'),
                         ukbb_coloc_cadd %>% select(VPRank) %>% mutate(Cat = 'CADD'))

ggplot(ukbb_coloc_types, aes(x=Cat,y=VPRank)) + geom_boxplot() + theme_bw() +
  ylab('Rank of variant p-value') + xlab('') +
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=16),
        legend.text=element_text(size=16),
        legend.title=element_blank()) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))



