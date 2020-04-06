library(data.table)
library(dplyr)
library(ggplot2)

data_dir = '/Users/nicoleferraro/durga_local/data/goats_data/v8_data/'

eqtl_coloc_results = fread(paste0(data_dir, 'coloc/coloc_sig_results_0.5.txt'))
eqtl_coloc_results$Trait = sapply(eqtl_coloc_results$filename, function(x) strsplit(x, '__PM__')[[1]][1])

eqtl_results = fread(paste0(data_dir, 'coloc/enloc_sig_results_0.5.txt'))
eqtl_results$Trait = sapply(eqtl_results$filename, function(x) strsplit(x, '_w_')[[1]][2])
eqtl_results$Trait = sapply(eqtl_results$Trait, function(x) strsplit(x, '_enloc_output.txt')[[1]][1])
sqtl_results = fread(paste0(data_dir, 'splicing/sqtl_enloc_rcp0.5_collapsed_bygene.txt'))
colnames(sqtl_results) = c('Chr', 'Start', 'End', 'rcp', 'Trait', 'Tissue', 'GeneChr', 'GeneStart', 'GeneEnd', 'Gene')
sqtl_results = sqtl_results %>% select(rcp, Trait, Gene)

## run below on cluster
all_watershed = fread(paste0('zcat ', data_dir, 'coloc/gtexV8.top.watershed.gam.per.type.per.variant.txt.gz'))
ws_data = data %>% filter(Model == 'Watershed') %>% group_by(Gene,Type) %>% top_n(1,ws_posterior) %>% ungroup()
gam_data = data %>% filter(Model == 'GAM') %>% group_by(Gene,Type) %>% top_n(1,gam_posterior) %>% ungroup()
##
ws_data_per_gene = fread(paste0(data_dir, 'coloc/gtexV8.top.watershed.gam.per.type.per.gene.txt'))
ws_data_per_gene$IsColoc = sapply(ws_data_per_gene$Gene, function(x) ifelse(x %in% coloc_genes, 1, 0))
ws_data_per_gene$IsColoc = factor(ws_data_per_gene$IsColoc, levels=c(0,1))

coloc_rr = data.frame(Risk = numeric(), Lower = numeric(), Upper = numeric(), Pval = numeric(), Type = character(), Model = character())
for (gtype in unique(ws_data_per_gene$Type)) {
  nv = length(unique(filter(ws_data_per_gene, Type == gtype, Model == 'Watershed', ws_posterior > 0.9)$VID))
  no_no = length(unique(filter(ws_data_per_gene, Type == gtype, Model == 'Watershed', ws_posterior < 0.9, IsColoc==0)$VID))
  no_yes = length(unique(filter(ws_data_per_gene, Type == gtype, Model == 'Watershed', ws_posterior < 0.9, IsColoc==1)$VID))
  yes_no = length(unique(filter(ws_data_per_gene, Type == gtype, Model == 'Watershed', ws_posterior > 0.9, IsColoc==0)$VID))
  yes_yes = length(unique(filter(ws_data_per_gene, Type == gtype, Model == 'Watershed', ws_posterior > 0.9, IsColoc==1)$VID))
  counttable = rbind(c(no_no, no_yes), c(yes_no, yes_yes))
  ws_rr = epitab(counttable,method="riskratio")$tab
  coloc_rr = rbind(coloc_rr, data.frame(Risk = ws_rr[2,5], Lower = ws_rr[2,6], Upper = ws_rr[2,7], Pval = ws_rr[2,8], Type = gtype, Model = 'Watershed'))
  
  gt_data = ws_data_per_gene %>% filter(Type == gtype, Model == 'GAM') %>% group_by(VID) %>% top_n(1,gam_posterior) %>% sample_n(1) %>% ungroup() %>% top_n(nv,gam_posterior)
  gt = min(gt_data$gam_posterior)
  no_no = length(unique(filter(ws_data_per_gene, Type == gtype, Model == 'GAM', ws_posterior < gt, IsColoc==0)$VID))
  no_yes = length(unique(filter(ws_data_per_gene, Type == gtype, Model == 'GAM', ws_posterior < gt, IsColoc==1)$VID))
  yes_no = length(unique(filter(ws_data_per_gene, Type == gtype, Model == 'GAM', ws_posterior >= gt, IsColoc==0)$VID))
  yes_yes = length(unique(filter(ws_data_per_gene, Type == gtype, Model == 'GAM', ws_posterior >= gt, IsColoc==1)$VID))
  counttable = rbind(c(no_no, no_yes), c(yes_no, yes_yes))
  gam_rr = epitab(counttable,method="riskratio")$tab
  coloc_rr = rbind(coloc_rr, data.frame(Risk = gam_rr[2,5], Lower = gam_rr[2,6], Upper = gam_rr[2,7], Pval = gam_rr[2,8], Type = gtype, Model = 'GAM'))
}

tcols = c('#7F5A83', '#BFCDE0', '#0D324D')
names(tcols) = c('ASE', 'TE', 'Splicing')
ggplot(coloc_rr, aes(x=Model,y=Risk,Group=Type)) + 
  geom_point(size=4,aes(color=Type),position=position_dodge(width=0.25)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width=0, position=position_dodge(width=0.25)) +
  geom_hline(yintercept=1,color='grey') + xlab('') + 
  scale_color_manual(values=tcols) +
  ylab('Relative risk of high posterior variant in coloc gene') +
  theme_bw() + theme(axis.text=element_text(size=16),
                     axis.title=element_text(size=18),
                     legend.title=element_blank(),
                     legend.position=c(0.1,0.1),
                     legend.text=element_text(size=16)) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


### Intersect ukbb stats
ukbb_sum_stats = fread(paste0(data_dir, 'coloc/gwas/gtex.coloc.traits.beta.gwas.overlap.txt'))
ukbb_sum_stats$Chr = sapply(ukbb_sum_stats$variant, function(x) strsplit(x, ':')[[1]][1])
ukbb_sum_stats$Pos = sapply(ukbb_sum_stats$variant, function(x) strsplit(x, ':')[[1]][2])
ukbb_sum_stats$VID = paste(paste0('chr',ukbb_sum_stats$Chr), ukbb_sum_stats$Pos, sep=':')

trait_table = fread(paste0(data_dir, 'coloc/gwas/ukbb.trait.table.txt'))
colnames(trait_table)[1] = 'Trait'
ukbb_sum_stats = inner_join(ukbb_sum_stats, trait_table, by='Trait')
ukbb_sum_stats = ukbb_sum_stats %>% group_by(Trait) %>% mutate(VRank = rank(abs(beta))) %>% ungroup()

### integrate watershed
exp_ws = fread(paste0(data_dir, 'coloc/gwas/expression.gam.watershed.ukbb.variants.txt'))
exp_ws$Gene = sapply(exp_ws$sample_names, function(x) strsplit(x, ':')[[1]][2])
sp_ws = fread(paste0(data_dir, 'coloc/gwas/splicing.gam.watershed.ukbb.variants.txt'))
sp_ws$Gene = sapply(sp_ws$sample_names, function(x) strsplit(x, ':')[[1]][2])
ase_ws = fread(paste0(data_dir, 'coloc/gwas/ase.gam.watershed.ukbb.variants.txt'))
ase_ws$Gene = sapply(ase_ws$sample_names, function(x) strsplit(x, ':')[[1]][2])

exp_ws = inner_join(exp_ws, ukbb_sum_stats %>% select(minor_AF,low_confidence_variant,beta,pval,VRank,Trait,Description,VID), by='VID')
sp_ws = inner_join(sp_ws, ukbb_sum_stats %>% select(minor_AF,low_confidence_variant,beta,pval,VRank,Trait,Description,VID), by='VID')
ase_ws = inner_join(ase_ws, ukbb_sum_stats %>% select(minor_AF,low_confidence_variant,beta,pval,VRank,Trait,Description,VID), by='VID')

ws_all = rbind(exp_ws %>% mutate(outlier_pvalue = total_expression_outlier_pvalue, gam_posterior = total_expression_gam_posterior, ws_posterior = total_expression_watershed_posterior) %>% 
                 select(-total_expression_outlier_pvalue, -total_expression_gam_posterior, -total_expression_watershed_posterior) %>% mutate(Type = 'TE'),
               sp_ws %>% mutate(outlier_pvalue = splicing_outlier_pvalue, gam_posterior = splicing_gam_posterior, ws_posterior = splicing_watershed_posterior) %>% 
                 select(-splicing_outlier_pvalue, -splicing_gam_posterior, -splicing_watershed_posterior) %>% mutate(Type = 'Splicing'),
               ase_ws %>% mutate(outlier_pvalue = ase_outlier_pvalue, gam_posterior = ase_gam_posterior, ws_posterior = ase_watershed_posterior) %>% 
                 select(-ase_outlier_pvalue, -ase_gam_posterior, -ase_watershed_posterior) %>% mutate(Type = 'ASE'))

rm(exp_ws, sp_ws, ase_ws)

ukbb.coloc.map = fread(paste0(data_dir, 'coloc/gwas/ukbb.traits.coloc.enloc.sqtl.map.txt'))
ws_coloc = filter(ws_all, Trait %in% ukbb.coloc.map$Trait)
ws_coloc = inner_join(ws_coloc, ukbb.coloc.map, by=c('Trait', 'Description'))

eqtl_results = eqtl_results %>% select(rcp, Trait, gene_id)
colnames(eqtl_results) = c('ColocProb', 'Enloc_name', 'Gene')
eqtl_coloc_results = eqtl_coloc_results %>% select(p4,Trait,gene)
colnames(eqtl_coloc_results) = c('ColocProb', 'Coloc_name', 'Gene')
colnames(sqtl_results) = c('ColocProb', 'Sqtl_name', 'Gene')

ws_eqtl_coloc = inner_join(ws_coloc, eqtl_coloc_results, by=c('Coloc_name', 'Gene'))
ws_eqtl_enloc = inner_join(ws_coloc, eqtl_results, by=c('Enloc_name', 'Gene'))
ws_sqtl = inner_join(ws_coloc, sqtl_results, by=c('Sqtl_name', 'Gene'))

ws_colocalized_regions = rbind(ws_eqtl_coloc %>% select(-sample_names,-Coloc_name,-Enloc_name,-Sqtl_name),
                               ws_eqtl_enloc %>% select(-sample_names,-Coloc_name,-Enloc_name,-Sqtl_name),
                               ws_sqtl %>% select(-sample_names,-Coloc_name,-Enloc_name,-Sqtl_name))

cadd_scores = fread(paste0(data_dir, 'coloc/gwas/gtexV8.ukbb.rare.cadd.scores.txt'))
ws_colocalized_regions = inner_join(ws_colocalized_regions, cadd_scores %>% select(RawScore,VID),by='VID')

ws_sum_per_variant = ws_colocalized_regions %>% group_by(VID) %>%
  mutate(MaxWS = max(ws_posterior, na.rm=T)) %>%
  mutate(MaxGAM = max(gam_posterior, na.rm=T)) %>%
  mutate(MaxCADD = max(RawScore, na.rm=T)) %>% sample_n(1) %>% ungroup()


ws_all = inner_join(ws_all, cadd_scores %>% select(RawScore,VID),by='VID')

high_ws = filter(ws_colocalized_regions, ws_posterior > 0.9) %>%
  group_by(VID,Trait) %>% sample_n(1) %>% 
  ungroup() %>% mutate(Model = 'Watershed')
## to match high ws > 0.9 for all types, gam > 0.17 and CADD > 6
high_gam = filter(ws_colocalized_regions, gam_posterior > 0.17) %>%
  group_by(VID,Trait) %>% sample_n(1) %>% 
  ungroup() %>% mutate(Model = 'GAM')
high_cadd = filter(ws_colocalized_regions, RawScore > 6) %>%
  group_by(VID,Trait) %>% sample_n(1) %>% 
  ungroup() %>% mutate(Model = 'CADD')

ws_top_sum = melt(rbind(high_ws, high_gam, high_cadd) %>% select(VID,Gene,VRank,Description,Model), id.vars=c('VID', 'Gene', 'Description', 'Model'))

ggplot(ws_top_sum, aes(x=Model,y=value)) + geom_boxplot() +
  theme_bw() + ylab('Effect size ranks') + xlab('') +
  theme(axis.title = element_text(size=16),
        axis.text = element_text(size=14)) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

## compare best GAM, WS, coloc
best_coloc = filter(ws_colocalized_regions, !is.na(outlier_pvalue)) %>%
  group_by(Gene) %>% top_n(1,VRank) %>% sample_n(1) %>% ungroup()
best_ws = filter(ws_all, !is.na(outlier_pvalue), ws_posterior > 0.5, Trait %in% best_coloc$Trait) %>%
  group_by(Gene) %>% top_n(1,VRank) %>% sample_n(1) %>% ungroup()
best_gam_vids = filter(ws_all, !is.na(outlier_pvalue)) %>% group_by(VID) %>%
  top_n(1,gam_posterior) %>% sample_n(1) %>% ungroup() %>%
  top_n(nrow(best_ws),gam_posterior) %>% select(VID)
best_gam = filter(ws_all, VID %in% best_gam_vids$VID, Trait %in% best_coloc$Trait) %>% group_by(Gene) %>%
  top_n(1,VRank) %>% sample_n(1) %>% ungroup()

all_ranks = rbind(best_ws %>% select(VRank) %>% mutate(Model = 'Watershed'),
                  best_gam %>% select(VRank) %>% mutate(Model = 'GAM'),
                  best_coloc %>% select(VRank) %>% mutate(Model = 'Colocalized'))

ggplot(all_ranks, aes(x=Model,y=VRank)) + geom_boxplot() + theme_bw()




