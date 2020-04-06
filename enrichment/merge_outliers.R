library(data.table)
library(dplyr)
library(ggplot2)

data_dir = '/users/nferraro/data/goats_data/v8_data/'

all_medz_data = fread(paste0(data_dir, 'gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.inds.genes.medz.txt'))
all_ase_data = fread(paste0(data_dir, 'ASE/combined.uncorrected.unfiltered.ad.median.scores.txt'))
splice_data = fread(paste0(data_dir, 'splicing/cross_tissue_covariate_method_none_no_global_outliers_ea_only_emperical_pvalue_gene_level.txt'))
splice_data = melt(splice_data)
all_medz_data$GeneID = sapply(all_medz_data$Gene, function(x) strsplit(x,'[.]')[[1]][1])
splice_data$GeneID = sapply(splice_data$CLUSTER_ID, function(x) strsplit(x,'[.]')[[1]][1])
colnames(all_ase_data)[1:3] = c('GeneID','Ind','MedDOT')
colnames(splice_data)[2:4] = c('Ind','MedP','GeneID')
splice_data$LogP = -log10(splice_data$MedP)
all_ase_data$LogDOT = -log10(all_ase_data$MedDOT)
splice_data$LogP = sapply(splice_data$LogP, function(x)
  ifelse(is.infinite(x), 6.3, x))
all_ase_data$LogDOT = sapply(all_ase_data$LogDOT, function(x)
  ifelse(is.infinite(x), 316.6, x))
medz_ase_data = merge(all_medz_data,all_ase_data,by=c('Ind','GeneID'))
medz_ase_splice_data = merge(medz_ase_data,splice_data,by=c('Ind','GeneID'))

write.table(medz_ase_splice_data,file=paste0(data_dir, 'gtexV8.medz.splice.ase.merged.txt'),sep='\t',quote=F,row.names=F)

medz_ase_splice_data = fread(paste0(data_dir, 'gtexV8.medz.splice.ase.merged.txt')) %>% select(Ind,GeneID,MedZ,LogDOT,LogP)
medz_splice = filter(medz_ase_splice_data,!is.na(LogDOT), !is.na(LogP))
rinds = sample(1:nrow(medz_splice),100000)
ggplot(medz_splice[rinds,],aes(x=LogP,y=LogDOT)) + geom_point() + theme_bw() +
  xlab('-log10(median splicing pval)') + ylab('-log10(median ASE pval)') +
  theme(axis.title = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.text.x = element_text(size=20),
        title = element_text(size=20)) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))




