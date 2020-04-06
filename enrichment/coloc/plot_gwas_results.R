library(data.table)
library(dplyr)
library(ggplot2)

gdir = '/srv/scratch/restricted/GOATs/features_v8/byGene/10kb_genebody/'
data_dir = '/users/nferraro/data/goats_data/v8_data/'

gtex_rare_vars = fread(paste0(gdir, 'all_rare_variants_SNPs.txt'),header=F) %>% select(-V5)
colnames(gtex_rare_vars) = c('Ind','Gene','Chr','Pos')
medz_outliers = fread(paste0(data_dir, 'gtexV8.outliers.v8ciseQTLs.removed.medz.Z2.txt'))
medz_outliers = merge(medz_outliers, gtex_rare_vars, by=c('Gene','Ind'))
splice_data = fread(paste0(data_dir, 'splicing/cross_tissue_covariate_method_none_no_global_outliers_ea_only_emperical_pvalue_gene_level.txt'))
splice_outliers = melt(splice_data) %>% filter(value < 2*pnorm(-abs(2)))
colnames(splice_outliers)[1:2] = c('Gene','Ind')
splice_outliers = merge(splice_outliers, gtex_rare_vars, by=c('Gene','Ind'))
rm(gtex_rare_vars)
write.table(medz_outliers, file=paste0(data_dir, 'coloc/gtexV8.outliers.v8ciseQTLs.removed.medz.withvars.txt'), sep='\t',quote=F,row.names=F)
write.table(splice_outliers, file=paste0(data_dir, 'coloc/gtexV8.outliers.v8ciseQTLs.removed.splicing.withvars.txt'), sep='\t',quote=F,row.names=F)

# medz_outliers = fread(paste0(data_dir, 'gtexV8.outliers.v8ciseQTLs.globalOutliers.removed.medz.withvars.txt'))
# wgs_inds = fread(paste0(data_dir, 'gtexv8_wgs_inds.txt'),header=F)
# globals = fread(paste0(data_dir, 'gtexV8_global_outliers_medz3_iqr.txt'),header=F)
# wgs_inds = filter(wgs_inds, !(V1 %in% globals$V1))
# 
# alltraits = fread(paste0(data_dir, 'coloc/gwas/sig.gtexV8.medz.enloc.overlap.txt'))
# alltraits$Chr = sapply(alltraits$variant, function(x) paste0('chr',strsplit(x,':')[[1]][1]))
# alltraits$Pos = sapply(alltraits$variant, function(x) strsplit(x,':')[[1]][2])
# 
# medz_outliers$Pos = as.character(medz_outliers$Pos)
# outlier_traits = merge(medz_outliers, alltraits, by=c('Chr','Pos'),all.x=T)
# outliers_sig_traits = filter(outlier_traits, !is.na(variant))
# 
# outs_per_ind = medz_outliers %>% filter(Ind %in% wgs_inds$V1) %>% group_by(Ind) %>% mutate(NumPerInd = length(unique(Gene))) %>% sample_n(1)
# rvs_per_ind = medz_outliers %>% filter(Ind %in% wgs_inds$V1) %>%
#   group_by(Ind,Gene) %>% mutate(HasVar = ifelse(any(!is.na(Chr)), 1, 0)) %>% sample_n(1) %>% ungroup() %>% group_by(Ind) %>% mutate(NumPerInd = sum(HasVar)) %>% sample_n(1) %>% ungroup()
# rvs_gwas_per_ind = outlier_traits %>% filter(Ind %in% wgs_inds$V1) %>%
#   group_by(Ind,Gene) %>% mutate(HasVar = ifelse(any(!is.na(variant)), 1, 0)) %>% sample_n(1) %>% ungroup() %>% group_by(Ind) %>% mutate(NumPerInd = sum(HasVar)) %>% sample_n(1) %>% ungroup()
# 
# coloc_results = fread(paste0(data_dir, 'coloc/enloc_sig_results_0.5.txt'))
# colnames(coloc_results)[1] = 'Gene'
# medz_coloc = merge(medz_outliers, coloc_results %>% select(Gene,rcp,lead_snp,lead_snp_rcp), by='Gene',all.x=T)
# coloc_per_ind = medz_coloc %>% filter(Ind %in% wgs_inds$V1) %>%
#   group_by(Ind,Gene) %>% mutate(HasColoc = ifelse(any(!is.na(lead_snp)), 1, 0)) %>% sample_n(1) %>% ungroup() %>% group_by(Ind) %>% mutate(NumPerInd = sum(HasColoc)) %>% sample_n(1) %>% ungroup()
# 
# all_per_ind = rbind(as.data.frame(outs_per_ind %>% select(Ind,NumPerInd) %>% mutate(Type = 'Outliers')), 
#                     as.data.frame(rvs_per_ind %>% select(Ind, NumPerInd) %>% mutate(Type = 'Rare_variants')),
#                     as.data.frame(coloc_per_ind %>% select(Ind, NumPerInd) %>% mutate(Type = 'Colocalization')),
#                     as.data.frame(rvs_gwas_per_ind %>% select(Ind, NumPerInd) %>% mutate(Type = 'Rare_trait_associated')))
# all_per_ind$Type = factor(all_per_ind$Type, levels=c('Outliers', 'Rare_variants', 'Colocalization', 'Rare_trait_associated'))
# ggplot(all_per_ind, aes(x=Type,y=NumPerInd)) + geom_boxplot() + theme_bw() +
#   xlab('') + ylab('Number per individual') +
#   theme(axis.title = element_text(size=20),
#         axis.text = element_text(size=20)) +
#   theme(panel.border = element_blank(),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"))
# 
# 
# trait_table = fread(paste0(data_dir, 'coloc/gwas/phenotype_codes.txt'))
# colnames(trait_table)[1] = 'Trait'
# alltraits = merge(alltraits, trait_table, by='Trait',all.x=T)
# outliers_sig_traits = merge(outliers_sig_traits, trait_table, by=c('Trait'),all.x=T)
# outliers_sig_traits = outliers_sig_traits %>% group_by(Trait,Chr,Pos,Gene,Ind) %>%
#   sample_n(1) %>% ungroup()
# 
# platelet_data = fread(paste0(data_dir, 'coloc/gwas/30100_irnt.gtexV8.medz.enloc.overlap.txt'))
# ggplot(platelet_data, aes(minor_AF)) + geom_histogram() + theme_bw() +
#   xlab('UKBB MAF') +
#   theme(axis.title = element_text(size=20),
#         axis.text = element_text(size=20)) +
#   theme(panel.border = element_blank(),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"))
# 
# sig_traits = filter(alltraits, pval < 0.05)
# sig_traits = merge(sig_traits, trait_table, by='Trait')
# sig_traits = merge(sig_traits, enloc_vars %>% select(variant, trait, tissue, prob, Ind, Gene), by='variant')
# 
# hypo_data = fread(paste0(data_dir, 'gwas/20002_1226.gwas.imputed_v3.both_sexes.tsv'))
# gstart = 637293 - 50000
# gend = 640706 + 50000
# hypo_data$Chr = sapply(hypo_data$variant, function(x) as.numeric(strsplit(x, ':')[[1]][1]))
# hypo_data$Pos = sapply(hypo_data$variant, function(x) as.numeric(strsplit(x, ':')[[1]][2]))
# sig_hypo = filter(hypo_data, Chr == 11, Pos > gstart & Pos < gend)
# vid = '11:648256:A:G'
# eid = '11:599313:T:C'
# sig_hypo$VOI = rep(0,nrow(sig_hypo))
# sig_hypo$VOI[which(sig_hypo$variant == vid)] = 1
# sig_hypo$VOI[which(sig_hypo$variant == eid)] = 2
# sig_hypo$VOI = factor(sig_hypo$VOI, levels=c(2,1,0))
# ggplot(sig_hypo, aes(x=minor_AF,y=beta,Group=VOI)) + geom_point(aes(color=VOI,size=VOI)) +
#   theme_bw() + scale_color_manual(values=c('green','red', 'black')) +
#   scale_size_manual(values=c(5,5,2)) +
#   guides(size=F,color=F) +
#   theme(axis.title = element_text(size=20),
#         axis.text = element_text(size=20)) +
#   theme(panel.border = element_blank(),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"))
# 
