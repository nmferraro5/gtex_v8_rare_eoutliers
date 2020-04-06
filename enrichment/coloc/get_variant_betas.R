library(data.table)
library(dplyr)
library(ggplot2)

gdir = '/srv/scratch/restricted/GOATs/features_v8/byGene/10kb_genebody/'
data_dir = '/users/nferraro/data/goats_data/v8_data/'

# gtex_rare_vars = fread(paste0(gdir, 'all_rare_variants_SNPs.txt'),header=F) %>% select(-V5)
# colnames(gtex_rare_vars) = c('Ind','Gene','Chr','Pos')
# medz_data = fread(paste0(data_dir, 'gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.medz.filtered.txt'))
# medz_outliers = merge(filter(medz_data, abs(MedZ) > 2), gtex_rare_vars, by=c('Gene','Ind')) %>%
#   mutate(VID = paste(Chr,Pos,sep=':'))
# medz_controls = filter(medz_data, abs(MedZ) < 1, Gene %in% medz_outliers$Gene)
# non_outlier_vars = merge(medz_controls, gtex_rare_vars, by=c('Gene','Ind')) %>%
#   mutate(VID = paste(Chr,Pos,sep=':')) %>% filter(!(VID %in% medz_outliers$VID))
# 
# splice_data = fread(paste0(data_dir, 'splicing/cross_tissue_covariate_method_none_no_global_outliers_ea_only_emperical_pvalue_gene_level.txt'))
# splice_melted = melt(splice_data)
# rm(splice_data)
# colnames(splice_melted)[1:2] = c('Gene','Ind')
# splice_outliers = filter(splice_melted, value < 2*pnorm(-abs(2)))
# splice_outliers = merge(splice_outliers, gtex_rare_vars, by=c('Gene','Ind')) %>%
#   mutate(VID = paste(Chr,Pos,sep=':'))
# splice_controls = filter(splice_melted, value > 2*pnorm(-abs(1)), Gene %in% splice_outliers$Gene)
# non_splice_outlier_vars = merge(splice_controls, gtex_rare_vars, by=c('Gene','Ind')) %>%
#   mutate(VID = paste(Chr,Pos,sep=':')) %>% filter(!(VID %in% splice_outliers$VID))
# rm(gtex_rare_vars)
# write.table(medz_outliers, file=paste0(data_dir, 'coloc/gtexV8.outliers.v8ciseQTLs.removed.medz.withvars.txt'), sep='\t',quote=F,row.names=F)
# write.table(splice_outliers, file=paste0(data_dir, 'coloc/gtexV8.outliers.v8ciseQTLs.removed.splicing.withvars.txt'), sep='\t',quote=F,row.names=F)
# write.table(non_outlier_vars, file=paste0(data_dir, 'coloc/gtexV8.controls.v8ciseQTLs.removed.medz.unmatchedvars.txt'), sep='\t',quote=F,row.names=F)
# write.table(non_splice_outlier_vars, file=paste0(data_dir, 'coloc/gtexV8.controls.v8ciseQTLs.removed.splicing.unmatchedvars.txt'), sep='\t',quote=F,row.names=F)

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

medz_vars_outliers = fread(paste0(data_dir, 'coloc/gtexV8.outliers.v8ciseQTLs.removed.medz.withvars.txt'),data.table=F) %>%
  filter(Gene %in% coloc_results$Gene)
splice_vars_outliers = fread(paste0(data_dir, 'coloc/gtexV8.outliers.v8ciseQTLs.removed.splicing.withvars.txt'),data.table=F) %>%
  filter(Gene %in% sqtl_data$Gene)

## compare beta distribution
medz_unmatched = fread(paste0(data_dir, 'coloc/gtexV8.controls.v8ciseQTLs.removed.medz.unmatchedvars.txt')) %>%
  filter(Gene %in% coloc_results$Gene)
splice_unmatched = fread(paste0(data_dir, 'coloc/gtexV8.controls.v8ciseQTLs.removed.splicing.unmatchedvars.txt')) %>%
  filter(Gene %in% sqtl_data$Gene)

top_traits = fread(paste0(data_dir, 'coloc/gwas/top.irnt.beta.gtexV8.medz.gwas.overlap.txt'))
top_traits$Chr = sapply(top_traits$variant, function(x) paste0('chr',strsplit(x,':')[[1]][1]))
top_traits$Pos = sapply(top_traits$variant, function(x) strsplit(x,':')[[1]][2])
top_traits$VID = paste(top_traits$Chr, top_traits$Pos, sep=':')

medz_vars_traits = merge(medz_vars_outliers, top_traits, by='VID') %>%
  mutate(Value = MedZ) %>% select(VID,Gene,Ind,Value,beta,pval,Trait,minor_AF) %>%
  mutate(Type = 'TotalExpression', Cat = 'Outlier')
splice_vars_traits = merge(splice_vars_outliers, top_traits, by='VID') %>%
  mutate(Value = value) %>% select(VID,Gene,Ind,Value,beta,pval,Trait,minor_AF) %>%
  mutate(Type = 'Splicing', Cat = 'Outlier')
medz_un_traits = merge(medz_unmatched, top_traits, by='VID') %>%
  mutate(Value = MedZ) %>% select(VID,Gene,Ind,Value,beta,pval,Trait,minor_AF) %>%
  mutate(Type = 'TotalExpression', Cat = 'Non-outlier')
splice_un_traits = merge(splice_unmatched, top_traits, by='VID') %>%
  mutate(Value = value) %>% select(VID,Gene,Ind,Value,beta,pval,Trait,minor_AF) %>%
  mutate(Type = 'Splicing', Cat = 'Non-outlier')

all_betas = rbind(medz_vars_traits, splice_vars_traits, medz_un_traits, splice_un_traits)
write.table(all_betas, file=paste0(data_dir, 'coloc/gwas/medz.splicing.top.irnt.nomsig.beta.coloc.genes.txt'),sep='\t',quote=F,row.names=F)

all_betas = fread(paste0(data_dir, 'coloc/gwas/medz.splicing.top.irnt.beta.coloc.genes.txt'))
all_betas = all_betas %>% group_by(VID,beta,Type) %>% sample_n(1) %>% ungroup()
ke_genes = filter(all_betas, Type == 'TotalExpression', Cat == 'Non-outlier')$Gene
ks_genes = filter(all_betas, Type == 'Splicing', Cat == 'Non-outlier')$Gene
te_betas = filter(all_betas, Type == 'TotalExpression', Gene %in% ke_genes)
splice_betas = filter(all_betas, Type == 'Splicing', Gene %in% ks_genes)
te_outlier_controls = filter(te_betas, Cat == 'Outlier')
for (gene in unique(te_outlier_controls$Gene)) {
  out_vars = filter(te_betas, Cat == 'Outlier',Gene == gene)$VID
  for (vid in out_vars) {
    maf = unique(filter(te_betas, Cat == 'Outlier',VID == vid)$minor_AF)
    controls = filter(te_betas, Cat == 'Non-outlier', Gene == gene) %>%
      mutate(DAF = abs(maf - minor_AF)) %>% top_n(1,-DAF) %>% select(-DAF)
    te_outlier_controls = rbind(te_outlier_controls, controls)
  }
}

sp_outlier_controls = filter(splice_betas, Cat == 'Outlier')
for (gene in unique(sp_outlier_controls$Gene)) {
  out_vars = filter(splice_betas, Cat == 'Outlier',Gene == gene)$VID
  for (vid in out_vars) {
    maf = unique(filter(splice_betas, Cat == 'Outlier',VID == vid)$minor_AF)
    controls = filter(splice_betas, Cat == 'Non-outlier', Gene == gene) %>%
      mutate(DAF = abs(maf - minor_AF)) %>% top_n(1,-DAF) %>% select(-DAF)
    sp_outlier_controls = rbind(sp_outlier_controls, controls)
  }
}

plot_betas = rbind(te_outlier_controls, sp_outlier_controls)
plot_betas_Z2 = plot_betas %>% mutate(Z=2)
plot_betas_Z3 = plot_betas %>% mutate(Z=3)
plot_betas = rbind(plot_betas_Z2,plot_betas_Z3)
plot_betas$Z = factor(plot_betas$Z, levels=c(2,3))
ggplot(plot_betas, aes(x=Cat,y=abs(beta),Group=Cat)) + geom_boxplot() + theme_bw() +
  xlab('') + facet_grid(.~Type) +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.text=element_text(size=18),
        legend.title=element_blank(),
        strip.text.x = element_text(size=16))

pheno_table = fread(paste0(data_dir, 'coloc/gwas/phenotype_codes.txt'))
colnames(pheno_table)[1] = 'Trait'
pheno_table = pheno_table %>% group_by(Trait,Description) %>% sample_n(1) %>% ungroup()
all_betas = merge(all_betas, pheno_table, by='Trait',all.x=T)

sig_traits = fread(paste0(data_dir, 'coloc/gwas/sig.gtexV8.medz.enloc.overlap.txt'))
sig_traits$Chr = sapply(sig_traits$variant, function(x) paste0('chr',strsplit(x,':')[[1]][1]))
sig_traits$Pos = sapply(sig_traits$variant, function(x) strsplit(x,':')[[1]][2])
sig_traits$VID = paste(sig_traits$Chr, sig_traits$Pos, sep=':')

medz_sig = merge(medz_vars_outliers, sig_traits, by='VID') %>%
  filter(Gene %in% coloc_results$Gene)
splice_sig = merge(splice_vars_outliers, sig_traits, by='VID') %>%
  filter(Gene %in% sqtl_data$Gene)

medz_sig = merge(medz_sig, pheno_table, by='Trait') %>%
  group_by(VID,Gene,Ind,Description) %>% top_n(1,abs(beta)) %>% sample_n(1) %>% ungroup() %>%
  select(VID,Gene,Ind,MedZ,minor_AF,beta,pval,Trait,Description)
splice_sig = merge(splice_sig, pheno_table, by='Trait') %>%
  group_by(VID,Gene,Ind,Description) %>% top_n(1,abs(beta)) %>% sample_n(1) %>% ungroup() %>%
  select(VID,Gene,Ind,value,minor_AF,beta,pval,Trait,Description)

medz_sig = merge(medz_sig, coloc_results %>% select(Gene,tissue,rcp), by='Gene')
splice_sig = merge(splice_sig, sqtl_data %>% select(Gene,Trait,rcp), by='Gene')


### compare betas with top outlier+control per trait

medz_summary = fread(paste0(data_dir, 'coloc/gwas/top.trait.beta.gtexV8.medz.gwas.variant.overlap.filtered.summary.txt'),data.table=F)

## bmi example
bmi_data = fread(paste0(data_dir, 'coloc/gwas/21001_raw.gwas.imputed_v3.both_sexes_chr11.tsv'))
testing = colnames(bmi_data)
bmi_data = bmi_data %>% select(-pval)
colnames(bmi_data) = testing[-4]
bmi_data$Chr = sapply(bmi_data$variant, function(x) strsplit(x, ':')[[1]][1])
bmi_data = filter(bmi_data, Chr == '11')
bmi_data$Pos = sapply(bmi_data$variant, function(x) strsplit(x, ':')[[1]][2])
bmi_data$VID = paste(bmi_data$Chr, bmi_data$Pos, sep=':')
gstart = 832884 - 10000
gend = 839720 + 10000
bmi_data$Pos = as.numeric(bmi_data$Pos)
bmi_data_filtered = filter(bmi_data, Pos > gstart & Pos < gend)
vid = '11:825341'
eqtl = '11:833668'
bmi_data_filtered$Type = sapply(bmi_data_filtered$VID, function(x)
  ifelse(x == vid, 'Variant', 
         ifelse(x == eqtl, 'eQTL', 'other')))

bmi_data_filtered$Type = factor(bmi_data_filtered$Type, levels=c('other','eQTL', 'Variant'))
ggplot(bmi_data_filtered, aes(x=minor_AF,y=beta)) +
  geom_point(aes(color=Type,size=Type)) 

#### trait by trait
bmi_genes = unique(filter(coloc_results, tissue == 'BMI')$Gene)
height_genes = unique(filter(coloc_results, tissue == 'Standing' | tissue == 'Height')$Gene)
bmi_splice_genes = unique(filter(sqtl_data, Trait == 'UKB_21001_Body_mass_index_BMI')$Gene)
height_splice_genes = unique(filter(sqtl_data, Trait == 'UKB_50_Standing_height')$Gene)
  
bmi_outliers = filter(medz_vars_outliers, Gene %in% bmi_genes)
bmi_controls = filter(medz_unmatched, Gene %in% bmi_genes)

height_outliers = filter(medz_vars_outliers, Gene %in% height_genes)
height_controls = filter(medz_unmatched, Gene %in% height_genes)

bmi_splice_outliers = filter(splice_vars_outliers, Gene %in% bmi_splice_genes)
bmi_splice_controls = filter(splice_unmatched, Gene %in% bmi_splice_genes)
height_splice_outliers = filter(splice_vars_outliers, Gene %in% height_splice_genes)
height_splice_controls = filter(splice_unmatched, Gene %in% height_splice_genes)

bmi_rare_vars = fread(paste0(data_dir, 'coloc/gwas/21001_irnt.gtexV8.medz.enloc.overlap.txt'))
bmi_rare_vars$Chr = sapply(bmi_rare_vars$variant, function(x) strsplit(x, ':')[[1]][1])
bmi_rare_vars$Pos = sapply(bmi_rare_vars$variant, function(x) strsplit(x, ':')[[1]][2])
bmi_rare_vars$VID = paste(paste0('chr',bmi_rare_vars$Chr), bmi_rare_vars$Pos, sep=':')

bmi_medz_outliers = merge(splice_vars_outliers, bmi_rare_vars, by='VID')
bmi_medz_controls = merge(splice_unmatched, bmi_rare_vars, by='VID')

bmi_medz_both = rbind(bmi_medz_outliers %>% mutate(Y='outlier'),bmi_medz_controls %>% mutate(Y='control'))
bmi_medz_both = bmi_medz_both %>% 
  filter(low_confidence_variant == 'FALSE') %>%
  group_by(variant,Y) %>% sample_n(1) %>% ungroup() %>%
  arrange(by=abs(beta))

bmi_medz_both$variant = factor(bmi_medz_both$variant, levels=bmi_medz_both$variant)
ggplot(bmi_medz_both, aes(x=Y,y=abs(beta))) + geom_boxplot() + theme_bw() +
  ggtitle('BMI effect size of splice outlier + non-outlier rare variants') +
  annotate("text",x=1.5,y=0.04,label='p = 0.61',cex=5) +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=18),
        legend.text=element_text(size=18),
        title=element_text(size=18),
        legend.title=element_blank()) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggplot(bmi_medz_both, aes(x=variant,y=abs(beta),Group=Y)) + 
  geom_bar(stat='identity',aes(fill=Y)) +
  theme_bw() + theme(axis.text.x = element_blank()) +
  ggtitle('Effect size of variants nearby genes co-localizing with BMI') +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.text=element_text(size=18),
        title=element_text(size=20),
        legend.title=element_blank()) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

bmi_sp_outliers = merge(height_splice_outliers, bmi_rare_vars, by='VID')
bmi_sp_controls = merge(height_splice_controls, bmi_rare_vars, by='VID')

bmi_splice_both = rbind(bmi_sp_outliers %>% mutate(Y='outlier'),bmi_sp_controls %>% mutate(Y='control'))

bmi_splice_both = bmi_splice_both %>% 
  filter(low_confidence_variant == 'FALSE') %>%
  filter(minor_AF > 0.0016) %>% filter(minor_AF < 0.0068) %>%
  group_by(variant,Y) %>% sample_n(1) %>% ungroup() %>%
  arrange(by=abs(beta))

bmi_splice_both$variant = factor(bmi_splice_both$variant, levels=bmi_splice_both$variant)
ggplot(bmi_splice_both, aes(x=variant,y=abs(beta),Group=Y)) + 
  geom_bar(stat='identity',aes(fill=Y)) +
  theme_bw() + theme(axis.text.x = element_blank()) +
  ggtitle('Effect size of variants nearby genes co-localizing with BMI') +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.text=element_text(size=18),
        title=element_text(size=20),
        legend.title=element_blank()) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))




