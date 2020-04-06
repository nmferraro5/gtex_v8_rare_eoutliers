library(data.table)
library(dplyr)
library(ggplot2)

# data_dir = '/srv/scratch/nferraro/gtex_temp/'
# 
# gtex_gwas = fread(paste0(data_dir, 'gtex.traits.beta.gwas.overlap.txt'))
# gtex_gwas$Chr = sapply(gtex_gwas$variant, function(x) strsplit(x, ':')[[1]][1])
# gtex_gwas$Pos = sapply(gtex_gwas$variant, function(x) strsplit(x, ':')[[1]][2])
# gtex_gwas$VID = paste(paste0('chr',gtex_gwas$Chr), gtex_gwas$Pos, sep=':')
# trait_table = fread(paste0(data_dir, 'ukbb.trait.table.txt'))
# colnames(trait_table)[1] = 'Trait'
# gtex_gwas = inner_join(gtex_gwas, trait_table, by='Trait')
# gtex_gwas = gtex_gwas %>% group_by(VID,Trait) %>% sample_n(1) %>% ungroup()
# 
# gene_pos = fread(paste0(data_dir, 'gencode.v26.GRCh38.genes_padded10kb.bed'))
# colnames(gene_pos) = c('Chr', 'Start', 'End', 'Gene')
# 
# all_outlier_data = fread(paste0(data_dir, 'gtex.trait.tissue.outliers.Z3.txt'))
data_dir = '/users/nferraro/data/goats_data/v8_data/'
trait_map = fread(paste0(data_dir, 'coloc/gwas/gtex.gwas.tissue.specific.traits.txt'))

get_ts_data <- function(i) {
  gtrait = trait_map$Description[i]
  print(gtrait)
  trait_data = filter(gtex_gwas, Description == gtrait)
  gtissue = trait_map$Tissue[i]
  tissue_data = filter(all_outlier_data, Tissue == gtissue)
  trait_data$Chr = paste0('chr', trait_data$Chr)
  trait_data = inner_join(trait_data, gene_pos, by='Chr') 
  
  trait_data$Pos = as.numeric(trait_data$Pos)
  kinds = which(trait_data$Pos > trait_data$Start & trait_data$Pos < trait_data$End)
  trait_data = trait_data[kinds,]
  outlier_genes = unique(filter(tissue_data, abs(Zscore) > 3)$Gene)
  trait_data$IsOutlier = sapply(trait_data$Gene, function(x)
    ifelse(x %in% outlier_genes, 1, 0))
  trait_data = trait_data %>% group_by(VID) %>% top_n(1,IsOutlier) %>% 
    sample_n(1) %>% ungroup() 
  trait_data$Tissue = gtissue
  return(trait_data)
}

# all_tissue_data = do.call(rbind, lapply(1:nrow(trait_map), get_ts_data))
# write.table(all_tissue_data, file=paste0(data_dir, 'gtex.ukbb.tissue.specific.outlier.trait.betas.txt'), sep='\t', quote=F, row.names=F)

all_tissue_data = fread(paste0(data_dir, 'coloc/gwas/gtex.ukbb.tissue.specific.all.types.outlier.trait.betas.txt'))
all_tissue_data$IsOutlier = factor(all_tissue_data$IsOutlier, levels=c(0,1))
colnames(all_tissue_data)[7] = 'Code'
trait_map$Name[4] = 'Hypothyroidism'
all_tissue_data = inner_join(all_tissue_data, trait_map, by='Code')
all_tissue_data$IsOutlier = sapply(all_tissue_data$IsOutlier, function(x)
  ifelse(x == 1, 'Outlier gene', 'Non-outlier gene'))
tplot = ggplot(all_tissue_data, aes(x=Name, y=abs(beta))) + geom_boxplot(aes(fill=IsOutlier)) +
  theme_bw() + scale_y_log10() + xlab('') +
  scale_fill_manual(values=c('grey', 'darkcyan')) +
  annotate("text", x=1:13, y=0.05, label='****',cex=4) +
  theme(axis.text.x=element_text(size=12, hjust=1, angle=25),
        axis.text.y=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text=element_text(size=8),
        legend.title=element_blank(),
        legend.key.height=unit(0.08,'in'),
        legend.position = c(0.95, 0.95)) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

cor_outliers = fread(paste0(data_dir, 'tissue_specific_outliers/gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.subset.Z2.specific.outliersOnly.jun1.txt'))
cor_outliers$Tissue = sapply(cor_outliers$Specific.Group, function(x) ifelse(x == 'Aorta', 'Artery_Aorta',
                                                                             ifelse(x == 'Adipose', 'Adipose_Subcutaneous', x)))

all_tissue_data = all_tissue_data %>% filter(low_confidence_variant == 'FALSE') %>%
  select(VID,minor_AF,AC,beta,Code,Gene,IsOutlier,Tissue.x,Name)
colnames(all_tissue_data)[8] = c('Tissue')
cor_outliers$TID = paste(cor_outliers$Gene, cor_outliers$Tissue, sep=':')
all_tissue_data$TID = paste(all_tissue_data$Gene, all_tissue_data$Tissue, sep=':')
all_tissue_cor_data = all_tissue_data %>% filter(Tissue %in% cor_outliers$Tissue) %>%
  mutate(IsTSOutlier = ifelse(TID %in% cor_outliers$TID, 'TS outlier', 'Non-outlier'))
tplot = ggplot(all_tissue_cor_data, aes(x=Name, y=abs(beta))) + geom_boxplot(aes(fill=IsTSOutlier)) +
  theme_bw() + scale_y_log10() + xlab('') +
  scale_fill_manual(values=c('grey', 'darkcyan')) +
  theme(axis.text.x=element_text(size=12, hjust=1, angle=25),
        axis.text.y=element_text(size=12),
        axis.title=element_text(size=14),
        legend.text=element_text(size=8),
        legend.title=element_blank(),
        legend.key.height=unit(0.08,'in'),
        legend.position = c(0.95, 0.95)) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

for (gname in unique(all_tissue_cor_data$Name)) {
  wtest = wilcox.test(abs(filter(all_tissue_cor_data, Name == gname, IsTSOutlier == 'TS outlier')$beta),
              abs(filter(all_tissue_cor_data, Name == gname, IsTSOutlier != 'TS outlier')$beta),alternative='g')
  print(gname)
  print(wtest$p.value)
}

# for (gname in unique(all_tissue_data$Name)) {
#   print(gname)
#   outlier_data = filter(all_tissue_data, Name == gname, IsOutlier == 'Outlier gene')
#   con_data = filter(all_tissue_data, Name == gname, IsOutlier != 'Outlier gene')
#   print(wilcox.test(abs(outlier_data$beta), abs(con_data$beta), alternative='g')$p.value)
# }


median_tpm = fread('/users/nferraro/data/goats_data/v8_data/All.tissues.median.tpm.txt')
median_brain = filter(median_tpm, grepl('Brain', Tissue))
median_brain$Tissue = 'Brain'
median_brain = median_brain %>% group_by(Gene) %>% mutate(MedBrain = median(MedianTPM)) %>% sample_n(1) %>%
  ungroup() 
median_other = filter(median_tpm, Tissue %in% all_tissue_data$Tissue.x)
median_brain = median_brain %>% select(Gene,MedBrain,Tissue)
colnames(median_brain)[2] = 'MedianTPM'
median_tpm = rbind(median_other, median_brain)

outlier_data = filter(all_tissue_data, IsOutlier == 'Outlier gene') %>% 
  select(minor_AF,beta,VID,Description,Gene,IsOutlier,Tissue)
nonoutlier_data = filter(all_tissue_data, IsOutlier == 'Non-outlier gene') %>% 
  select(minor_AF,beta,VID,Description,Gene,IsOutlier,Tissue)
colnames(outlier_data)[7] = 'Tissue'
colnames(nonoutlier_data)[7] = 'Tissue'

outlier_data_tpm = inner_join(outlier_data, median_tpm, by=c('Tissue','Gene'))
nonoutlier_data_tpm = inner_join(nonoutlier_data, median_tpm, by=c('Tissue','Gene'))

outlier_genes = outlier_data_tpm %>% group_by(Tissue,Description,Gene) %>% sample_n(1) %>% ungroup()
nonoutlier_genes = nonoutlier_data_tpm %>% group_by(Tissue,Description,Gene) %>% sample_n(1) %>% ungroup()

both_data = inner_join(outlier_genes %>% select(Gene,Description,Tissue,MedianTPM), 
                       nonoutlier_genes %>% select(Gene,Description,Tissue,MedianTPM), by=c('Description','Tissue'))

both_data$TPM_Dif = abs(both_data$MedianTPM.x - both_data$MedianTPM.y)
both_data = both_data %>% group_by(Gene.x,Tissue,Name) %>% top_n(1,-TPM_Dif) %>% ungroup() %>%
  filter(TPM_Dif < 0.01)

both_data$OID = paste(both_data$Gene.x, both_data$Tissue, both_data$Name, sep='_')
both_data$NID = paste(both_data$Gene.y, both_data$Tissue, both_data$Name, sep='_')
all_tissue_data$UID = paste(all_tissue_data$Gene, all_tissue_data$Tissue.x, all_tissue_data$Name, sep='_')
tissue_outliers = filter(all_tissue_data, IsOutlier == 'Outlier gene', UID %in% both_data$OID)
tissue_nonoutliers = filter(all_tissue_data, IsOutlier == 'Non-outlier gene', UID %in% both_data$NID)
matched_tissue_data = rbind(tissue_outliers, tissue_nonoutliers)

matched_tissue_data = matched_tissue_data %>% group_by(Name) %>%
  mutate(VRank = rank(abs(beta))) %>% ungroup()

best_matched_tissue_data = matched_tissue_data %>% filter(!is.na(beta)) %>% group_by(Name,Gene) %>%
  top_n(1,VRank) %>% ungroup()
tremove = c('BMI', 'Body fat percentage', 'Neuroticism', 'Migraine self-reported', 'Schizophrenia', 'MS', 'Epilepsy diagnosis')
ggplot(best_matched_tissue_data %>% filter(!(Name %in% tremove)), aes(x=Name, y=abs(beta))) + geom_boxplot(aes(fill=IsOutlier)) +
  theme_bw() + xlab('') + scale_y_log10() +
  scale_fill_manual(values=c('grey', 'darkcyan')) +
  annotate("text", x=1:6, y=0.0075, label='****', cex=4) +
  theme(axis.text.x=element_text(size=13, hjust=1, angle=25),
        axis.text.y=element_text(size=13),
        axis.title=element_text(size=15),
        legend.text=element_text(size=12),
        legend.title=element_blank(),
        legend.key.height=unit(0.1,'in')) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggplot(best_matched_tissue_data, aes(x=IsOutlier,y=VRank)) + geom_violin(aes(fill=IsOutlier)) + 
  geom_boxplot(width=0.4) + 
  theme_bw() + xlab('') + ylab('Variant effect size rank') +
  scale_fill_manual(values=c('grey', 'darkcyan')) + guides(fill=F) +
  annotate("text", x=1.5, y=15000, label='****', cex=7) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16)) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

for (gname in unique(matched_tissue_data$Name)) {
  print(gname)
  outlier_data = filter(best_matched_tissue_data, Name == gname, IsOutlier == 'Outlier gene')
  con_data = filter(best_matched_tissue_data, Name == gname, IsOutlier != 'Outlier gene')
  #print(wilcox.test(abs(outlier_data$beta), abs(con_data$beta), alternative='g')$p.value)
  print(wilcox.test(outlier_data$VRank, con_data$VRank, alternative='g')$p.value)
}

## integrate coloc data - come back to this
gdir = '/Users/nicoleferraro/durga_local/data/goats_data/v8_data/'
eqtl_coloc_results = fread(paste0(gdir, 'coloc/coloc_sig_results_0.5.txt'))
eqtl_coloc_results$Trait = sapply(eqtl_coloc_results$filename, function(x) strsplit(x, '__PM__')[[1]][1])

eqtl_results = fread(paste0(gdir, 'coloc/enloc_sig_results_0.5.txt'))
eqtl_results$Trait = sapply(eqtl_results$filename, function(x) strsplit(x, '_w_')[[1]][2])
eqtl_results$Trait = sapply(eqtl_results$Trait, function(x) strsplit(x, '_enloc_output.txt')[[1]][1])
sqtl_results = fread(paste0(gdir, 'splicing/sqtl_enloc_rcp0.5_collapsed_bygene.txt'))
colnames(sqtl_results) = c('Chr', 'Start', 'End', 'rcp', 'Trait', 'Tissue', 'GeneChr', 'GeneStart', 'GeneEnd', 'Gene')
sqtl_results = sqtl_results %>% select(rcp, Trait, Gene)

ukbb.coloc.map = fread(paste0(gdir, 'coloc/gwas/only.ukbb.traits.coloc.enloc.sqtl.map.txt'))
eqtl_results = eqtl_results %>% select(rcp, Trait, gene_id)
colnames(eqtl_results) = c('ColocProb', 'Enloc_name', 'Gene')
eqtl_coloc_results = eqtl_coloc_results %>% select(p4,Trait,gene)
colnames(eqtl_coloc_results) = c('ColocProb', 'Coloc_name', 'Gene')
colnames(sqtl_results) = c('ColocProb', 'Sqtl_name', 'Gene')

eqtl_results = inner_join(eqtl_results, ukbb.coloc.map, by='Enloc_name')
sqtl_results = inner_join(sqtl_results, ukbb.coloc.map, by='Sqtl_name')
eqtl_coloc_results = inner_join(eqtl_coloc_results, ukbb.coloc.map, by='Coloc_name')

coloc_all = rbind(eqtl_results %>% select(Gene,Description),
                  sqtl_results %>% select(Gene,Description),
                  eqtl_coloc_results %>% select(Gene,Description))

best_matched_coloc_data = best_matched_tissue_data %>% select(IsOutlier,VRank,VID,Gene,Name,Tissue.x,Tissue.y,Description.x)
colnames(best_matched_coloc_data)[8] = 'Description'
best_matched_coloc_data = inner_join(best_matched_coloc_data, coloc_all, by=c('Gene','Description'))

ggplot(best_matched_coloc_data, aes(x=IsOutlier,y=VRank)) + geom_violin(aes(fill=IsOutlier)) + 
  geom_boxplot(width=0.4) + 
  theme_bw() + xlab('') + ylab('Variant effect size rank') +
  scale_fill_manual(values=c('grey', 'darkcyan')) + guides(fill=F) +
  annotate("text", x=1.5, y=15000, label='****', cex=7) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16)) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))



