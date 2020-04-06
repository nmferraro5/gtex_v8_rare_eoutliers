library(data.table)
library(dplyr)
library(epitools)
library(ggplot2)

data_dir = '/users/nferraro/data/goats_data/v8_data/'
args = commandArgs(trailingOnly = T)
window = args[1]

## Define arguments

window_exp_data = fread(paste0(data_dir, 'enrichments/windows/ExpressionNeighboringOutliers_', window, '_Window.txt')) %>%
  mutate(DataType = 'TE')
window_ase_data = fread(paste0(data_dir, 'enrichments/windows/ASENeighboringOutliers_', window, '_Window.txt')) %>%
  mutate(DataType = 'ASE')

gencode = fread(paste0(data_dir, 'gencode.v26.GRCh38.genes.bed'))
gencode$Gene = sapply(gencode$V4, function(x) strsplit(x, '[.]')[[1]][1])
ase_genes = unique(c(window_ase_data$Gene1, window_ase_data$Gene2))
ase_genes = sapply(ase_genes, function(x) filter(gencode, Gene == x)$V4[1])
window_ase_data$Gene1 = sapply(window_ase_data$Gene1, function(x) filter(gencode, Gene == x)$V4[1])
window_ase_data$Gene2 = sapply(window_ase_data$Gene2, function(x) filter(gencode, Gene == x)$V4[1])

te.variants = fread(paste0(data_dir, 'enrichments/windows/TE.window.variant.cats.txt')) %>%
  mutate(DataType = 'TE')
ase.variants = fread(paste0(data_dir, 'enrichments/windows/ASE.window.variant.cats.txt')) %>%
  mutate(DataType = 'ASE')
all.variants = rbind(te.variants, ase.variants)
all.variants$Subject_ID = paste0('GTEX-', all.variants$indiv_id)

all_window_data = rbind(window_exp_data, window_ase_data)
colnames(all_window_data)[3] = 'gene_id'
all_window_data = merge(all_window_data, all.variants, by=c('gene_id', 'Subject_ID', 'DataType'), all.x=T)
colnames(all_window_data)[1] = 'Gene1'
colnames(all_window_data)[5] = 'gene_id'
all_window_data = merge(all_window_data, all.variants, by=c('gene_id', 'Subject_ID', 'DataType'), all.x=T)
colnames(all_window_data)[1] = 'Gene2'
all_window_data = all_window_data %>% select(Subject_ID, Gene1, Gene2, DataType, chr_numeric.x, pos.x, variant_cat.x, chr_numeric.y, pos.y, variant_cat.y)
variant_order = as.data.frame(cbind(1:14, c('DUP', 'CNV', 'DEL', 'BND', 'INV', 'TE', 'splice', 'frameshift', 'stop', 'TSS', 'conserved_noncoding', 'coding', 'other_noncoding', 'no_variant')))
rownames(variant_order) = variant_order[,2]
variant_order$V1 = as.numeric(as.character(variant_order$V1))
all_window_data$VCat1 = sapply(all_window_data$variant_cat.x, function(x) variant_order[x,1])
all_window_data$VCat2 = sapply(all_window_data$variant_cat.y, function(x) variant_order[x,1])

window_summary = all_window_data %>% group_by(Gene1,Gene2,Subject_ID,DataType) %>%
  mutate(BestCat1 = min(VCat1)) %>% mutate(BestCat2 = min(VCat2)) %>% top_n(1,-VCat1) %>%
  sample_n(1) %>% ungroup() %>% filter(!is.na(BestCat1)) %>% filter(!is.na(BestCat2))

window_summary = window_summary %>% mutate(Neither = ifelse(BestCat1 %in% c(13,14) & BestCat2 %in% c(13,14), 1, 0)) %>%
  mutate(JustOne = ifelse((BestCat1 %in% c(13,14) | BestCat2 %in% c(13,14)) & Neither == 0, 1, 0)) %>%
  group_by(Subject_ID,Gene1,Gene2,DataType) %>%
  mutate(BestOverall = min(BestCat1, BestCat2)) %>% sample_n(1) %>% ungroup()

window_cat1_summary = window_summary %>% filter(Neither == 0) %>% filter(JustOne == 0) %>%
  group_by(DataType) %>% mutate(N=n()) %>% ungroup() %>% group_by(DataType,BestCat1) %>%
  mutate(N1 = n()) %>% sample_n(1) %>% ungroup() 
window_cat2_summary = window_summary %>% filter(Neither == 0) %>% filter(JustOne == 0) %>%
  group_by(DataType) %>% mutate(N=n()) %>% ungroup() %>% group_by(DataType,BestCat2) %>%
  mutate(N2 = n()) %>% sample_n(1) %>% ungroup() 

## ASE outliers
ase.outliers = fread(paste0(data_dir, 'ASE/combined.ad.scores.in.MEDIAN_4_10_update.outliers.tsv'))
ase.outliers$Gene = sapply(ase.outliers$V1, function(x) filter(gencode, Gene == x)$V4[1])
ase.summary = filter(window_summary, DataType == 'ASE')
colnames(ase.outliers)[c(2,4)] = c('Subject_ID', 'Gene1')
ase.summary = inner_join(ase.summary, ase.outliers %>% select(-V1), by=c('Subject_ID', 'Gene1'))
colnames(ase.outliers)[c(2,4)] = c('Subject_ID', 'Gene2')
ase.summary = inner_join(ase.summary, ase.outliers %>% select(-V1), by=c('Subject_ID', 'Gene2'))

splice.summary = filter(window_summary, DataType == 'Splicing')
splice.summary$GenePair = paste(splice.summary$Gene1, splice.summary$Gene2, sep='_')
clusters = fread(paste0(data_dir, 'splicing/cluster_info.txt'))
genes = strsplit(clusters$genes, split = ",")
clusters_df = data.frame(cluster_id = rep(clusters$cluster_id, sapply(genes, length)), Gene = unlist(genes))
colnames(clusters_df)[2] = 'Gene1'
cluster.summary = inner_join(splice.summary, clusters_df, by='Gene1')
colnames(clusters_df)[2] = 'Gene2'
cluster.summary = inner_join(cluster.summary, clusters_df, by='Gene2')

# te.summary = filter(window_summary, DataType == 'TE')
# te.summary$GenePair = paste(te.summary$Gene1, te.summary$Gene2, sep='_')
# te.overlap = te.summary %>% filter(BestCat1 == BestCat2) %>% filter(Neither == 0) %>%
#   select(-Neither,-JustOne,-VCat1,-VCat2,-chr_numeric.x, -chr_numeric.y,-pos.x,-pos.y,-variant_cat.y,-BestCat1,-BestCat2)
# 
# all_window_data$GenePair = paste(all_window_data$Gene1, all_window_data$Gene2, sep='_')
# 
# te.overlap = inner_join(te.overlap, all_window_data, by=c('Subject_ID', 'GenePair', 'variant_cat.x', 'DataType')) %>%
#   filter(pos.x == pos.y) %>% group_by(Subject_ID, GenePair, pos.x) %>% sample_n(1) %>% ungroup()
# 
# te.var.sum = te.overlap %>% group_by(Subject_ID, GenePair) %>% top_n(1,-BestOverall) %>% sample_n(1) %>% ungroup()
# var.sum = as.data.frame(table(te.var.sum$variant_cat.x))
# var.sum$VarProp = var.sum$Freq / 74
# var.sum$Var1 = factor(var.sum$Var1, levels=c('CNV', 'DUP', 'DEL', 'TSS', 'conserved_noncoding'))
# ggplot(var.sum, aes(x=Var1, y=VarProp)) + geom_bar(stat='identity') + theme_bw() +
#   xlab('') + ylab('Proportion of outlier pairs with same variant') +
#   theme(panel.grid.minor = element_blank(),panel.background = element_blank(), 
#         panel.border = element_blank(), panel.grid.major = element_blank(),
#         axis.text=element_text(size=14),
#         axis.title=element_text(size=16),
#         axis.text.x=element_text(hjust=1,angle=45))

## calculate type enrichments
te.outliers = fread(paste0(data_dir, 'gtexV8.outliers.v8ciseQTLs.globalOutliers.removed.medz.txt'))
colnames(te.outliers)[c(1,2)] = c('Subject_ID', 'gene_id')

all.variants$VCat = sapply(all.variants$variant_cat, function(x) variant_order[x,1])
sum.variants = all.variants %>% group_by(gene_id,indiv_id,DataType) %>% top_n(1,-VCat) %>%
  sample_n(1) %>% ungroup()
sum.variants$UID = paste(sum.variants$Subject_ID, sum.variants$gene_id, sep='_')

get_cat_rr <- function(dtype, vid) {
  vdata = filter(window_summary, DataType == dtype, BestOverall == vid)
  vdata$UID1 = paste(vdata$Subject_ID, vdata$Gene1, sep='_')
  vdata$UID2 = paste(vdata$Subject_ID, vdata$Gene2, sep='_')
  
  type.variants = filter(sum.variants, DataType == dtype)
  kgenes = unique(c(filter(window_summary, DataType == dtype)$Gene1, filter(window_summary, DataType == dtype)$Gene2))
  controls = filter(type.variants, gene_id %in% kgenes, !(UID %in% vdata$UID1), !(UID %in% vdata$UID2))
  if (dtype == 'ASE') {
    colnames(ase.outliers)[4] = 'gene_id'
    controls = inner_join(controls, ase.outliers, by=c('Subject_ID', 'gene_id'))
  } else {
    controls = inner_join(controls, te.outliers, by=c('Subject_ID', 'gene_id'))
  }
  
  yes_yes = nrow(vdata)
  yes_no = nrow(filter(window_summary, DataType == dtype, BestOverall != vid))
  no_no = nrow(filter(controls, VCat != vid))
  no_yes = nrow(filter(controls, VCat == vid))
  counttable = rbind(c(no_no, no_yes), c(yes_no, yes_yes))
  cat.rr = epitab(counttable, method='riskratio')$tab
  return(data.frame(VCat = vid, DataType = dtype, Riskratio = cat.rr[2,5], Lower = cat.rr[2,6], Upper = cat.rr[2,7], Pval = cat.rr[2,8]))
}

all.ase.rrs = do.call(rbind, lapply(variant_order$V1, function(x) get_cat_rr('ASE', x)))
all.te.rrs = do.call(rbind, lapply(variant_order$V1, function(x) get_cat_rr('TE', x)))

all.ase.rrs$Cat = variant_order$V2
all.te.rrs$Cat = variant_order$V2

all.rrs = rbind(all.ase.rrs, all.te.rrs) %>% mutate(Window = window)
write.table(all.rrs, file=paste0(data_dir, 'enrichments/windows/ase.te.pair.enrichments.', window, '.txt'), sep='\t', quote=F, row.names=F)

all_files = dir(paste0(data_dir, 'enrichments/windows/'), 'ase.te.pair.enrichments.*.txt', full=T)
all.rrs = do.call(rbind,lapply(all_files, function(x) fread(x)))
all.rrs = filter(all.rrs, !(VCat %in% c(13,14))) %>% filter(!is.na(Upper)) %>% filter(!is.na(Lower))
all.rrs$Window = factor(all.rrs$Window, levels=c('100kbp', '500kbp', '1000kbp', '5000kbp'))


ggplot(all.rrs %>% filter(Window == '100kbp', DataType == 'TE', Cat %in% c('CNV', 'DUP', 'TSS')), aes(x=Cat, y=Riskratio)) + 
  geom_point(size=3) + xlab('') + ylab('Relative risk') +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width=0) +
  theme_bw() + scale_y_log10() + geom_hline(yintercept=1, linetype="dashed", colour="darkgrey")








