library(data.table)
library(dplyr)
library(ggplot2)

gdir = '/srv/scratch/restricted/GOATs/features_v8/byGene/10kb_genebody/'
data_dir = '/users/nferraro/data/goats_data/v8_data/'

gtex_rare_vars = fread(paste0(gdir, 'all_rare_variants_SNPs.txt'),header=F) %>% select(-V5)
colnames(gtex_rare_vars) = c('Ind','Gene','Chr','Pos')
medz_data = fread(paste0(data_dir, 'gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.medz.filtered.txt'))
medz_outliers = merge(filter(medz_data, abs(MedZ) > 2), gtex_rare_vars, by=c('Gene','Ind')) %>%
  mutate(VID = paste(Chr,Pos,sep=':'))
medz_controls = filter(medz_data, abs(MedZ) < 1, Gene %in% medz_outliers$Gene)
non_outlier_vars = merge(medz_controls, gtex_rare_vars, by=c('Gene','Ind')) %>%
  mutate(VID = paste(Chr,Pos,sep=':')) %>% filter(!(VID %in% medz_outliers$VID))

splice_data = fread(paste0(data_dir, 'splicing/cross_tissue_covariate_method_none_no_global_outliers_ea_only_emperical_pvalue_gene_level.txt'))
splice_melted = melt(splice_data)
rm(splice_data)
colnames(splice_melted)[1:2] = c('Gene','Ind')
splice_outliers = filter(splice_melted, value < 2*pnorm(-abs(2)))
splice_outliers = merge(splice_outliers, gtex_rare_vars, by=c('Gene','Ind')) %>%
  mutate(VID = paste(Chr,Pos,sep=':'))
splice_controls = filter(splice_melted, value > 2*pnorm(-abs(1)), Gene %in% splice_outliers$Gene)
non_splice_outlier_vars = merge(splice_controls, gtex_rare_vars, by=c('Gene','Ind')) %>%
  mutate(VID = paste(Chr,Pos,sep=':')) %>% filter(!(VID %in% splice_outliers$VID))
rm(gtex_rare_vars)
write.table(medz_outliers, file=paste0(data_dir, 'coloc/gtexV8.outliers.v8ciseQTLs.removed.medz.withvars.txt'), sep='\t',quote=F,row.names=F)
write.table(splice_outliers, file=paste0(data_dir, 'coloc/gtexV8.outliers.v8ciseQTLs.removed.splicing.withvars.txt'), sep='\t',quote=F,row.names=F)
write.table(non_outlier_vars, file=paste0(data_dir, 'coloc/gtexV8.controls.v8ciseQTLs.removed.medz.unmatchedvars.txt'), sep='\t',quote=F,row.names=F)
write.table(non_splice_outlier_vars, file=paste0(data_dir, 'coloc/gtexV8.controls.v8ciseQTLs.removed.splicing.unmatchedvars.txt'), sep='\t',quote=F,row.names=F)

### get number per individual
wgs_inds = fread(paste0(data_dir, 'gtexv8_wgs_inds.txt'),header=F)
globals = fread(paste0(data_dir, 'gtexV8_global_outliers_medz3_iqr.txt'),header=F)
wgs_inds = filter(wgs_inds, !(V1 %in% globals$V1))

medz_outliers = fread(paste0(data_dir, 'gtexV8.outliers.v8ciseQTLs.removed.medz.Z2.txt')) %>%
  filter(Ind %in% wgs_inds$V1)
splice_outliers = fread(paste0(data_dir, 'splicing/gtexV8.outliers.v8ciseQTLs.removed.splicing.Z2.txt')) %>%
  filter(Ind %in% wgs_inds$V1)
coloc_results = fread(paste0(data_dir, 'coloc/enloc_sig_results_0.5.txt'))
colnames(coloc_results)[1] = 'Gene'
sqtl_data = fread(paste0(data_dir, 'splicing/sqtl_enloc_rcp0.5_collapsed_bygene.txt'))
sqtl_data = sqtl_data %>% select(-V7,-V8,-V9)
colnames(sqtl_data) = c('Chr','Start','End','rcp','Trait','Tissue','Gene')
medz_outliers$Ind = factor(medz_outliers$Ind, levels=unique(medz_outliers$Ind))
medz_all_per_ind_2 = as.data.frame(table(medz_outliers$Ind))
medz_all_per_ind_3 = as.data.frame(table(filter(medz_outliers, abs(MedZ) > 3)$Ind))

medz_per_ind_2 = as.data.frame(table(filter(medz_outliers, Gene %in% coloc_results$Gene)$Ind))
medz_per_ind_3 = as.data.frame(table(filter(medz_outliers, Gene %in% coloc_results$Gene, abs(MedZ) > 3)$Ind))

splice_outliers$Ind = factor(splice_outliers$Ind, levels=unique(splice_outliers$Ind))
splice_per_ind_2 = as.data.frame(table(filter(splice_outliers, Gene %in% sqtl_data$Gene)$Ind))
splice_per_ind_3 = as.data.frame(table(filter(splice_outliers, Gene %in% sqtl_data$Gene, value < 0.0027)$Ind))
splice_all_per_ind_2 = as.data.frame(table(splice_outliers$Ind))
splice_all_per_ind_3 = as.data.frame(table(filter(splice_outliers, value < 0.0027)$Ind))

medz_vars_outliers = fread(paste0(data_dir, 'coloc/gtexV8.outliers.v8ciseQTLs.removed.medz.withvars.txt'),data.table=F) %>%
  filter(Gene %in% coloc_results$Gene)
splice_vars_outliers = fread(paste0(data_dir, 'coloc/gtexV8.outliers.v8ciseQTLs.removed.splicing.withvars.txt'),data.table=F) %>%
  filter(Gene %in% sqtl_data$Gene)
medz_vars_sum2 = medz_vars_outliers %>% group_by(Ind) %>% mutate(NRV = length(unique(Gene))) %>%
  sample_n(1) %>% select(Ind,NRV) %>% ungroup()
medz_vars_sum3 = medz_vars_outliers %>% filter(abs(MedZ) > 3) %>% group_by(Ind) %>%
  mutate(NRV = length(unique(Gene))) %>%
  sample_n(1) %>% select(Ind,NRV) %>% ungroup()
splice_vars_sum2 = splice_vars_outliers %>% group_by(Ind) %>% mutate(NRV = length(unique(Gene))) %>%
  sample_n(1) %>% select(Ind,NRV) %>% ungroup()
splice_vars_sum3 = splice_vars_outliers %>% filter(value < 0.0027) %>% group_by(Ind) %>%
  mutate(NRV = length(unique(Gene))) %>%
  sample_n(1) %>% select(Ind,NRV) %>% ungroup()

medz_all_vars_sum2 = medz_all_vars_outliers %>% group_by(Ind) %>% mutate(NRV = length(unique(Gene))) %>%
  sample_n(1) %>% select(Ind,NRV) %>% ungroup()
medz_all_vars_sum3 = medz_all_vars_outliers %>% filter(abs(MedZ) > 3) %>% group_by(Ind) %>%
  mutate(NRV = length(unique(Gene))) %>%
  sample_n(1) %>% select(Ind,NRV) %>% ungroup()
splice_all_vars_sum2 = splice_all_vars_outliers %>% group_by(Ind) %>% mutate(NRV = length(unique(Gene))) %>%
  sample_n(1) %>% select(Ind,NRV) %>% ungroup()
splice_all_vars_sum3 = splice_all_vars_outliers %>% filter(value < 0.0027) %>% group_by(Ind) %>%
  mutate(NRV = length(unique(Gene))) %>%
  sample_n(1) %>% select(Ind,NRV) %>% ungroup()

medz_per_ind_2 = merge(medz_per_ind_2 %>% mutate(Ind=Var1) %>% select(-Var1), medz_vars_sum2, by='Ind',all.x=T)
medz_per_ind_3 = merge(medz_per_ind_3 %>% mutate(Ind=Var1) %>% select(-Var1), medz_vars_sum3, by='Ind',all.x=T)
splice_per_ind_2 = merge(splice_per_ind_2 %>% mutate(Ind=Var1) %>% select(-Var1), splice_vars_sum2, by='Ind',all.x=T)
splice_per_ind_3 = merge(splice_per_ind_3 %>% mutate(Ind=Var1) %>% select(-Var1), splice_vars_sum3, by='Ind',all.x=T)

medz_all_per_ind_2 = merge(medz_all_per_ind_2 %>% mutate(Ind=Var1) %>% select(-Var1), medz_all_vars_sum2, by='Ind',all.x=T)
medz_all_per_ind_3 = merge(medz_all_per_ind_3 %>% mutate(Ind=Var1) %>% select(-Var1), medz_all_vars_sum3, by='Ind',all.x=T)
splice_all_per_ind_2 = merge(splice_all_per_ind_2 %>% mutate(Ind=Var1) %>% select(-Var1), splice_all_vars_sum2, by='Ind',all.x=T)
splice_all_per_ind_3 = merge(splice_all_per_ind_3 %>% mutate(Ind=Var1) %>% select(-Var1), splice_all_vars_sum3, by='Ind',all.x=T)


splice_per_ind = rbind(splice_all_per_ind_2 %>% mutate(Z=2,Type='Splicing',Genes='All'),splice_per_ind_2 %>% mutate(Z=2, Type = 'Splicing',Genes='Colocalized'), splice_all_per_ind_3 %>% mutate(Z=3,Type='Splicing',Genes='All'), splice_per_ind_3 %>% mutate(Z=3, Type = 'Splicing',Genes='Colocalized'))
medz_per_ind = rbind(medz_all_per_ind_2 %>% mutate(Z=2,Type='TotalExpression',Genes='All'),medz_per_ind_2 %>% mutate(Z=2, Type = 'TotalExpression',Genes='Colocalized'), medz_all_per_ind_3 %>% mutate(Z=3,Type='TotalExpression',Genes='All'), medz_per_ind_3 %>% mutate(Z=3, Type = 'TotalExpression',Genes='Colocalized'))
both_per_ind = melt(rbind(medz_per_ind, splice_per_ind), id.vars=c('Ind','Z','Type','Genes'))
both_per_ind$variable = sapply(both_per_ind$variable, function(x)
  ifelse(x == 'Freq', 'Outliers', 'Outliers+RV'))
both_per_ind$value = sapply(both_per_ind$value, function(x) ifelse(is.na(x), 0, x))
both_per_ind$Z = factor(both_per_ind$Z, levels=c(2,3))

both_per_ind$variable2 = paste(both_per_ind$Genes, both_per_ind$variable,sep='_')
both_per_ind$variable2 = factor(both_per_ind$variable2, levels=c('All_Outliers', 'All_Outliers+RV', 'Colocalized_Outliers', 'Colocalized_Outliers+RV'))
ggplot(both_per_ind, aes(x=variable2,y=value+1,Group=Z)) + geom_boxplot(aes(fill=Z)) + theme_bw() +
  xlab('') + ylab('Log(outliers/ind + 1)') + facet_grid(.~Type) +
  scale_y_log10() + scale_fill_manual(values=c('gray30','gray66')) +
  theme(axis.text=element_text(size=18),
        axis.text.x=element_text(size=16,hjust=1,angle=35),
        axis.title=element_text(size=20),
        legend.text=element_text(size=18),
        legend.title=element_text(size=18),
        strip.text.x = element_text(size=16))

## look at sig overlap
all_genes = fread(paste0(data_dir, 'gencode.v26.GRCh38.genes_genetypes_autosomal_PCandlinc_only.txt'),header=F)
em = length(unique(filter(coloc_results, Gene %in% all_genes$V1)$Gene))
en = nrow(filter(all_genes, !(V1 %in% coloc_results$Gene)))
ek2 = length(unique(medz_outliers$Gene))
ex2 = length(which(unique(medz_outliers$Gene) %in% coloc_results$Gene))
medz_hg2 = phyper(ex2,em,en,ek2,lower.tail=F)
ek3 = length(unique(filter(medz_outliers,abs(MedZ)>3)$Gene))
ex3 = length(which(unique(filter(medz_outliers,abs(MedZ)>3)$Gene) %in% coloc_results$Gene))
medz_hg3 = phyper(ex3,em,en,ek3,lower.tail=F)
m = length(unique(filter(sqtl_data, Gene %in% all_genes$V1)$Gene))
n = nrow(filter(all_genes, !(V1 %in% sqtl_data$Gene)))
k2 = length(unique(splice_outliers$Gene))
x2 = length(which(unique(splice_outliers$Gene) %in% sqtl_data$Gene))
splice_hg2 = phyper(x2,m,n,k2,lower.tail=F)
k3 = length(unique(filter(splice_outliers,value<0.0027)$Gene))
x3 = length(which(unique(filter(splice_outliers,value < 0.0027)$Gene) %in% sqtl_data$Gene))
splice_hg3 = phyper(x3,m,n,k3,lower.tail=F)

overlap_data = as.data.frame(rbind(c(em,ek2,ex2,2,'TotalExpression'),
                     c(em,ek3,ex3,3,'TotalExpression'),
                     c(m,k2,x2,2,'Splicing'),
                     c(m,k3,x3,3,'Splicing')))
colnames(overlap_data) = c('Colocalized_Genes','Outliers','Overlap','ZT','Type')
overlap_melt = melt(overlap_data,id.vars=c('ZT','Type'))
overlap_melt$value = as.numeric(as.character(overlap_melt$value))
overlap_melt$variable = factor(overlap_melt$variable, levels=c('Overlap','Outliers','Colocalized_Genes'))
ggplot(overlap_melt, aes(x=ZT,y=value,Group=variable)) +
  geom_bar(stat='identity',aes(fill=variable)) +
  facet_grid(.~Type) + theme_bw() +
  scale_fill_manual(values=c('#C7D3DD','#7D82B8', '#613F75')) +
  xlab('Outlier threshold') + ylab('Number of genes') +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.title=element_blank(),
        legend.text=element_text(size=16),
        strip.text.x=element_text(size=16))

## compare beta distribution
medz_unmatched = fread(paste0(data_dir, 'coloc/gtexV8.controls.v8ciseQTLs.removed.medz.unmatchedvars.colocgenes.txt'))
splice_unmatched = fread(paste0(data_dir, 'coloc/gtexV8.controls.v8ciseQTLs.removed.splice.unmatchedvars.colocgenes.txt'))

top_traits = fread(paste0(data_dir, 'coloc/gwas/top.irnt.nomsig.beta.gtexV8.medz.gwas.overlap.txt'))
top_traits$Chr = sapply(top_traits$variant, function(x) paste0('chr',strsplit(x,':')[[1]][1]))
top_traits$Pos = sapply(top_traits$variant, function(x) strsplit(x,':')[[1]][2])
top_traits$VID = paste(top_traits$Chr, top_traits$Pos, sep=':')

medz_vars_traits = merge(medz_vars_outliers, top_traits, by='VID')
splice_vars_traits = merge(splice_vars_outliers, top_traits, by='VID')
medz_un_traits = merge(medz_unmatched, top_traits, by='VID') 
splice_un_traits = merge(splice_unmatched, top_traits, by='VID')

outlier_top = 0
control_top = 0
for (cgene in unique(medz_vars_traits$Gene)) {
  outlier_beta = max(abs(filter(medz_vars_traits, Gene == cgene)$beta))
  control_beta = sample(1,abs(filter(medz_un_traits, Gene == cgene)$beta))
  if (outlier_beta > control_beta) { 
    outlier_top = outlier_top + 1
  } else {
    control_top = control_top + 1
  }
}

hypo_data = fread(paste0(data_dir, 'gwas/20002_1226.gwas.imputed_v3.both_sexes.tsv'))
gstart = 637293 - 50000
gend = 640706 + 50000
hypo_data$Chr = sapply(hypo_data$variant, function(x) as.numeric(strsplit(x, ':')[[1]][1]))
hypo_data$Pos = sapply(hypo_data$variant, function(x) as.numeric(strsplit(x, ':')[[1]][2]))
sig_hypo = filter(hypo_data, Chr == 11, Pos > gstart & Pos < gend)
vid = '11:648256:A:G'
eid = '11:599313:T:C'
sig_hypo$VOI = rep(0,nrow(sig_hypo))
sig_hypo$VOI[which(sig_hypo$variant == vid)] = 1
sig_hypo$VOI[which(sig_hypo$variant == eid)] = 2
sig_hypo$VOI = factor(sig_hypo$VOI, levels=c(2,1,0))
ggplot(sig_hypo, aes(x=minor_AF,y=beta,Group=VOI)) + geom_point(aes(color=VOI,size=VOI)) +
  theme_bw() + scale_color_manual(values=c('green','red', 'black')) +
  scale_size_manual(values=c(5,5,2)) +
  guides(size=F,color=F) +
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=20)) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))



