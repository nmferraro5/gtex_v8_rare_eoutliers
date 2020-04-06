library(data.table)
library(dplyr)
library(ggplot2)

data_dir = '/users/nferraro/data/goats_data/v8_data/'
gdir = '/srv/scratch/restricted/GOATs/features_v8/byGene/10kb_genebody/'

gtex_gwas = fread(paste0(data_dir, 'coloc/gwas/gtex.traits.beta.gwas.overlap.txt'))
gtex_gwas$Chr = sapply(gtex_gwas$variant, function(x) strsplit(x, ':')[[1]][1])
gtex_gwas$Pos = sapply(gtex_gwas$variant, function(x) strsplit(x, ':')[[1]][2])
gtex_gwas$VID = paste(paste0('chr',gtex_gwas$Chr), gtex_gwas$Pos, sep=':')

medz_outliers = fread(paste0(data_dir, 'gtexV8.outliers.v8ciseQTLs.removed.medz.Z2.txt')) %>%
  mutate(GID = paste(Ind,Gene,sep='_'))
medz_outliers_wrare = fread(paste0(data_dir, 'gtexV8.outliers.v8ciseQTLs.removed.medz.Z2.withrarevars.txt'))
medz_outliers_wrare = medz_outliers_wrare %>% mutate(VID = paste(Chr,Pos,sep=':'))
medz_outliers_ukbb_vars = filter(medz_outliers_wrare, VID %in% gtex_gwas$VID)

sp_outliers = fread(paste0(data_dir, 'splicing/gtexV8.outliers.v8ciseQTLs.removed.splicing.Z2.txt')) %>%
  mutate(GID = paste(Ind,Gene,sep='_'))
sp_outliers_wrare = fread(paste0(data_dir, 'splicing/gtexV8.outliers.v8ciseQTLs.removed.splicing.Z2.withrarevars.txt'))
sp_outliers_wrare = sp_outliers_wrare %>% mutate(VID = paste(Chr,Pos,sep=':'))
sp_outliers_ukbb_vars = filter(sp_outliers_wrare, VID %in% gtex_gwas$VID) %>%
  mutate(GID = paste(Ind,Gene,sep='_'))

ase_data = fread(paste0(data_dir, 'ASE/combined.ad.scores.in.MEDIAN_4_10_update.tsv'))
ase_outliers = melt(ase_data) %>% filter(value < 2*pnorm(-abs(2))) 
colnames(ase_outliers)[1] = 'Gene'
gene_names = fread(paste0(data_dir, 'gencode.v26.GRCh38.genes.bed')) %>%
  select(-V1,-V2,-V3) %>% mutate(Gene = V4)
gene_names$Gene = sapply(gene_names$Gene, function(x) strsplit(x, '[.]')[[1]][1])
ase_outliers = inner_join(ase_outliers, gene_names, by='Gene') %>%
  mutate(GID = paste(variable,V4, sep='_'))

ase_outliers_wrare = fread(paste0(data_dir, 'ASE/gtexV8.ase.outliers.Z2.withrarevars.txt'))
ase_outliers_wrare = ase_outliers_wrare %>% mutate(VID = paste(Chr,Pos,sep=':'))
ase_outliers_ukbb_vars = filter(ase_outliers_wrare, VID %in% gtex_gwas$VID) %>%
  mutate(GID = paste(Ind,Gene,sep='_'))

all_ukbb_rare = rbind(medz_outliers_ukbb_vars %>% select(MAF,VID) %>% mutate(DataType = 'TE'),
                      sp_outliers_ukbb_vars %>% select(MAF,VID) %>% mutate(DataType = 'Splicing'),
                      ase_outliers_ukbb_vars %>% select(MAF,VID) %>% mutate(DataType = 'ASE'))

all_ukbb_rare = all_ukbb_rare %>% group_by(VID,DataType) %>% sample_n(1) %>% ungroup()
all_ukbb_rare$MAF = as.numeric(all_ukbb_rare$MAF)
mafplot = ggplot(all_ukbb_rare, aes(MAF,Group=DataType)) + 
  geom_density(aes(fill=DataType),alpha=0.5) + theme_bw() +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.text=element_text(size=18),
        legend.title=element_blank())

## getting matched non-outlier variants
ase_non_outliers = melt(ase_data) %>% filter(value < 2*pnorm(-abs(1))) 
colnames(ase_non_outliers)[1] = 'Gene'
ase_non_outliers = inner_join(ase_non_outliers, gene_names, by='Gene') %>%
  mutate(GID = paste(variable, Gene, sep='_'))
print(ase_non_outliers[1,])
medz_non_outliers = fread(paste0(data_dir, 'gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.medz.txt')) %>%
  filter(abs(MedZ) < 1) %>% mutate(GID = paste(Ind,Gene,sep='_'))
sp_non_outliers = melt(fread(paste0(data_dir, 'splicing/cross_tissue_covariate_method_none_no_global_outliers_ea_only_emperical_pvalue_gene_level.txt'))) %>%
  filter(value < 2*pnorm(-abs(1))) %>% mutate(GID = paste(variable,CLUSTER_ID,sep='_'))

rare_snps = fread(paste0(gdir, 'all_rare_variants_SNPs_MAF.txt'),header=F,fill=T)
te_matched = filter(rare_snps, V2 %in% medz_outliers_ukbb_vars$Gene) %>%
  mutate(VID = paste(V3,V4,sep=':')) %>% filter(!(VID %in% medz_outliers_ukbb_vars$VID)) %>%
  mutate(GID = paste(V1,V2,sep='_')) %>% filter(GID %in% medz_non_outliers$GID) %>%
  mutate(DataType = 'TE')
sp_matched = filter(rare_snps, V2 %in% sp_outliers_ukbb_vars$Gene) %>%
  mutate(VID = paste(V3,V4,sep=':')) %>% filter(!(VID %in% sp_outliers_ukbb_vars$VID)) %>%
  mutate(GID = paste(V1,V2,sep='_')) %>% filter(GID %in% sp_non_outliers$GID) %>%
  mutate(DataType = 'Splicing')
ase_matched = filter(rare_snps, V2 %in% ase_outliers_ukbb_vars$Gene) %>%
  mutate(VID = paste(V3,V4,sep=':')) %>% filter(!(VID %in% ase_outliers_ukbb_vars$VID)) %>%
  mutate(GID = paste(V1,V2,sep='_')) %>% filter(GID %in% ase_non_outliers$GID) %>%
  mutate(DataType = 'ASE')

all_matched = rbind(te_matched, sp_matched, ase_matched) %>% filter(VID %in% gtex_gwas$VID)
write.table(all_matched, file=paste0(data_dir, 'gtexV8.all.data.types.inlier.ukbb.rare.snps.txt'), sep='\t',quote=F,row.names=F)

non_outlier_vars = fread(paste0(data_dir, 'gtexV8.all.data.types.non.outlier.ukbb.rare.snps.txt'))
non_outlier_vars = non_outlier_vars %>% filter(VID %in% gtex_gwas$VID)
colnames(non_outlier_vars)[1:5] = c('Ind','Gene','Chr','Pos','MAF')
all_outliers_ukbb_vars = rbind(medz_outliers_ukbb_vars %>% select(Ind,Gene,VID) %>% mutate(DataType='TE'),
                               sp_outliers_ukbb_vars %>% select(Ind,Gene,VID) %>% mutate(DataType='Splicing'),
                               ase_outliers_ukbb_vars %>% select(Ind,Gene,VID) %>% mutate(DataType='ASE'))
all_outliers_ukbb_vars = all_outliers_ukbb_vars %>% mutate(Y = 'outlier') %>%
  filter(Gene %in% non_outlier_vars$Gene) %>%
  group_by(Gene,VID,DataType) %>% sample_n(1) %>% ungroup()
non_outlier_vars = non_outlier_vars %>% mutate(Y = 'control') %>%
  filter(Gene %in% all_outliers_ukbb_vars$Gene) %>%
  group_by(Gene,VID,DataType) %>% sample_n(1) %>% ungroup()
both_vars = rbind(all_outliers_ukbb_vars %>% select(-Ind),
                  non_outlier_vars %>% select(-Ind,-Chr,-Pos,-GID,-MAF))

trait_table = fread(paste0(data_dir, 'coloc/gwas/ukbb.trait.table.txt'))
colnames(trait_table)[1] = 'Trait'
both_vars = inner_join(both_vars, gtex_gwas, by='VID')
#both_vars = inner_join(both_vars, trait_table, by='Trait')
both_vars_false = filter(both_vars, low_confidence_variant == 'FALSE')
outlier_vars_false = filter(both_vars_false, Y == 'outlier') %>%
  group_by(Gene,VID,DataType,Trait) %>% sample_n(1) %>% ungroup()
control_vars_false = filter(both_vars_false, Y == 'control') %>%
  group_by(Gene,VID,DataType,Trait) %>% sample_n(1) %>% ungroup()

outlier_vars_false = inner_join(outlier_vars_false, ws_all %>% select(Gene,VID,DataType,posterior,Trait), by=c('Gene','VID','DataType','Trait')) %>%
  filter(posterior > 0.05)
control_vars_false = inner_join(control_vars_false, ws_all %>% select(Gene,VID,DataType,posterior,Trait), by=c('Gene','VID','DataType','Trait'))

#coloc_data = fread(paste0(data_dir, 'coloc/enloc_sig_results_0.5.txt'))

joined_vars = inner_join(outlier_vars_false %>% filter(minor_AF < 0.1) %>% select(-Chr,-Pos,-low_confidence_variant),
                         control_vars_false %>% filter(minor_AF < 0.1) %>% select(-Chr,-Pos,-low_confidence_variant),
                         by=c('Gene','DataType','Trait','Description'))
joined_vars = joined_vars %>% select(-Trait,-posterior.x,-posterior.y)
joined_vars = joined_vars %>% group_by(Gene,DataType,Description) %>%
  mutate(DAF = abs(minor_AF.x - minor_AF.y)) %>%
  top_n(1,-DAF) %>% ungroup() %>% filter(DAF < 0.001) %>%
  select(Gene, Description, DataType, Y.x, beta.x, Y.y, beta.y, minor_AF.x, minor_AF.y, DAF)


## possible traits - Standing height, Non-cancer illness code, self-reported: hypothyroidism/myxoedema  
melt_vars = melt(joined_vars,id.vars=c('Gene','Description','DataType', 'DAF', 'Y.x', 'Y.y')) %>%
  group_by(Gene,Description,DataType) %>%
  mutate(HasBoth = ifelse('beta.x' %in% variable & 'beta.y' %in% variable, 1, 0)) %>%
  ungroup() %>% filter(HasBoth == 1) %>%
  group_by(Gene,DataType,variable) %>% top_n(1,abs(value)) %>% ungroup()

#group_by(Gene,DataType,variable,value) %>% sample_n(1) %>% ungroup()

#group_by(Gene,DataType,variable) %>% top_n(1,abs(value)) %>% ungroup()

beta_plot = filter(melt_vars, variable %in% c('beta.x','beta.y'))
beta_plot$Type = sapply(beta_plot$variable, function(x)
  ifelse(x == 'beta.y', 'Non-outlier', 'Outlier'))
ggplot(beta_plot, aes(x=DataType,y=abs(value),Group=Type)) +
  geom_boxplot(aes(fill=Type)) +
  theme_bw() + xlab('') + ylab('Max abs(beta) per gene') +
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=20),
        title = element_text(size=20),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        legend.key.size = unit(1.5, 'lines')) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


ggplot(beta_plot, aes(x=Type,y=abs(value),Group=Type)) +
  geom_boxplot(aes(fill=Type)) +
  theme_bw() + xlab('') + ylab('Max abs(beta) per gene') +
  guides(fill=F) +
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=20),
        title = element_text(size=20),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        legend.key.size = unit(1.5, 'lines')) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggplot(beta_plot, aes(value,Group=Type)) +
  geom_density(aes(fill=Type),alpha=0.5) +
  theme_bw() + xlab('Max Beta') + guides(fill=F) +
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=16),
        title = element_text(size=20),
        legend.text=element_text(size=20),
        legend.title=element_blank(),
        legend.key.size = unit(1.5, 'lines'),
        strip.text.y=element_text(size=16)) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  facet_grid(DataType~.)

te_count = c(0,0)
sp_count = c(0,0)
ase_count = c(0,0)
for (gtrait in unique(beta_plot$Description)) {
  te_data = merge(filter(beta_plot, Description == gtrait, DataType == 'TE',Type == 'Outlier') %>% select(Description,Gene,value),
                  filter(beta_plot, Description == gtrait, DataType == 'TE', Type != 'Outlier') %>% select(Description,Gene,value),
                  by='Description') %>% mutate(OutlierHigher = ifelse(abs(value.x) > abs(value.y), 1, 0))
  te_count = c(te_count[1] + sum(te_data$OutlierHigher), te_count[2] + (nrow(te_data) - sum(te_data$OutlierHigher)))
  sp_data = merge(filter(beta_plot, Description == gtrait, DataType == 'Splicing',Type == 'Outlier') %>% select(Description,Gene,value),
                  filter(beta_plot, Description == gtrait, DataType == 'Splicing', Type != 'Outlier') %>% select(Description,Gene,value),
                  by='Description') %>% mutate(OutlierHigher = ifelse(abs(value.x) > abs(value.y), 1, 0))
  sp_count = c(sp_count[1] + sum(sp_data$OutlierHigher), sp_count[2] + (nrow(sp_data) - sum(sp_data$OutlierHigher)))
  ase_data = merge(filter(beta_plot, Description == gtrait, DataType == 'ASE',Type == 'Outlier') %>% select(Description,Gene,value),
                  filter(beta_plot, Description == gtrait, DataType == 'ASE', Type != 'Outlier') %>% select(Description,Gene,value),
                  by='Description') %>% mutate(OutlierHigher = ifelse(abs(value.x) > abs(value.y), 1, 0))
  ase_count = c(ase_count[1] + sum(ase_data$OutlierHigher), ase_count[2] + (nrow(ase_data) - sum(ase_data$OutlierHigher)))
}

maf_plot = filter(melt_vars, variable %in% c('minor_AF.x','minor_AF.y'))
maf_plot$Type = sapply(maf_plot$variable, function(x)
  ifelse(x == 'minor_AF.y', 'Non-outlier', 'Outlier'))
ggplot(maf_plot, aes(x=DataType,y=value,Group=Type)) +
  geom_boxplot(aes(fill=Type)) +
  theme_bw() + xlab('') + ylab('UKBB minor_AF') +
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=20),
        title = element_text(size=20),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        legend.key.size = unit(1.5, 'lines')) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))



