#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(ggplot2)
require(RColorBrewer)

data_dir = '/users/nferraro/data/goats_data/v8_data/'
# 
# knn_outliers = fread(paste0(data_dir, 'gtexV8.outliers.v8ciseQTLs.globalOutliers.removed.subset.knn.txt'),data.table=F)
# knn_outliers = knn_outliers %>% top_n(3281,-FDR)
# medz_outliers = fread(paste0(data_dir, 'gtexV8.outliers.v8ciseQTLs.globalOutliers.removed.subset.medz.txt'),data.table=F)
# 
#enloc_results = fread(paste0(data_dir, 'coloc/enloc_results.txt'),data.table=F)
# enloc_results = fread(paste0(data_dir, 'coloc/coloc_results_p4.txt'),data.table=F)
# colnames(enloc_results)[1] = 'Gene'
# medz_enloc = merge(medz_outliers,enloc_results,by='Gene',all.x=T)
# knn_enloc = merge(knn_outliers,enloc_results,by='Gene',all.x=T)
# rm(enloc_results)
#medz_enloc_top = medz_enloc %>% group_by(Gene,Ind) %>% top_n(1,rcp) %>% ungroup() %>% select(Gene,Ind,trait,tissue,rcp)
#knn_enloc_top = knn_enloc %>% group_by(Gene,Ind) %>% top_n(1,rcp) %>% ungroup() %>% select(Gene,Ind,trait,tissue,rcp)
# medz_enloc_top = medz_enloc %>% group_by(Gene,Ind) %>% top_n(1,p4) %>% ungroup() %>% select(Gene,Ind,trait,tissue,p4)
# knn_enloc_top = knn_enloc %>% group_by(Gene,Ind) %>% top_n(1,p4) %>% ungroup() %>% select(Gene,Ind,trait,tissue,p4)
# 
# write.table(medz_enloc_top,file=paste0(data_dir, 'gtexV8.outliers.v8ciseQTLs.globalOutliers.removed.medz.coloc.txt'),sep='\t',quote=F,row.names=F)
# write.table(knn_enloc_top,file=paste0(data_dir, 'gtexV8.outliers.v8ciseQTLs.globalOutliers.removed.knn.coloc.txt'),sep='\t',quote=F,row.names=F)
# 
#medz_enloc = fread(paste0(data_dir, 'gtexV8.outliers.v8ciseQTLs.globalOutliers.removed.medz.enloc.nopos.txt'),data.table=F)
#print(agenes)
#medz_enloc_top = medz_enloc %>% group_by(Gene,Ind) %>% mutate(BF = agenes*rcp) %>%
#    filter(BF < 0.01) %>% ungroup()
#print(nrow(medz_enloc_top))
#medz_enloc_top = medz_enloc_top %>% select(-N,-Df,-Y,-MedZ,-lead_snp,-lead_snp_rcp)
#write.table(medz_enloc_top,file=paste0(data_dir, 'gtexV8.outliers.v8ciseQTLs.globalOutliers.removed.medz.top.enloc.nopos.txt'),sep='\t',quote=F,row.names=F)
#
#enloc_sum = medz_enloc %>% filter(rcp>0.5) %>% mutate(UCID = paste(tissue,trait,sep='_')) %>% group_by(Gene,Ind) %>% mutate(NCs = length(unique(UCID))) %>% mutate(NTiss = length(unique(tissue))) %>% mutate(NTraits = length(unique(trait))) %>% top_n(1,rcp) %>% ungroup()
#write.table(enloc_sum, file=paste0(data_dir, 'gtexV8.outliers.v8ciseQTLs.globalOutliers.removed.medz.top.enloc.summary.txt'),sep='\t',quote=F,row.names=F)
#
#enloc_knn_sum = knn_enloc %>% filter(rcp>0.5) %>% mutate(UCID = paste(tissue,trait,sep='_')) %>% group_by(Gene,Ind) %>% mutate(NCs = length(unique(UCID))) %>% mutate(NTiss = length(unique(tissue))) %>% mutate(NTraits = length(unique(trait))) %>% top_n(1,rcp) %>% ungroup()
#write.table(enloc_knn_sum, file=paste0(data_dir, 'gtexV8.outliers.v8ciseQTLs.globalOutliers.removed.knn.top.enloc.summary.txt'),sep='\t',quote=F,row.names=F)
#
#data_dir = '/Users/nicoleferraro/durga_local/data/goats_data/v8_data/'
# enloc_sum = fread(paste0(data_dir, 'gtexV8.outliers.v8ciseQTLs.globalOutliers.removed.medz.coloc.txt'))
# enloc_knn_sum = fread(paste0(data_dir, 'gtexV8.outliers.v8ciseQTLs.globalOutliers.removed.knn.coloc.txt'))
# 
# # change p4 to rcp for enloc
# enloc_sum_top = enloc_sum %>% mutate(UID=paste(Gene,Ind,sep='_')) %>%
#  group_by(UID) %>% top_n(1,rcp) %>% sample_n(1) %>% ungroup() %>%
#  mutate(Method='MEDZ')
# enloc_knn_sum_top = enloc_knn_sum %>% mutate(UID=paste(Gene,Ind,sep='_')) %>%
#  group_by(UID) %>% top_n(1,rcp) %>% sample_n(1) %>% ungroup() %>%
#  mutate(Method='COR')
# 
# meth.colors = brewer.pal(9, 'YlGnBu')[3:9]
# names(meth.colors) = c('EM','COR','MEAN','PMD','SOFT','STFZ','MEDZ')
# both_rcp = rbind(enloc_sum_top,enloc_knn_sum_top)
# ggplot(both_rcp,aes(rcp,Group=Method)) + geom_histogram(aes(fill=Method),alpha=0.7) +
#  theme_bw() + scale_fill_manual(values=meth.colors) +
#  xlab('Colocalization probability') +
#  theme(axis.text=element_text(size=20),
#        axis.title=element_text(size=20),
#        legend.text=element_text(size=20),
#        legend.title=element_blank(),
#        panel.border = element_blank(),
#        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#        panel.background = element_blank(), axis.line = element_line(colour = "black"))
# 
# medz_top = filter(enloc_sum,rcp>0.5)
# knn_top = filter(enloc_knn_sum,rcp>0.5)

# get trait counts
enloc_sum = fread(paste0(data_dir, 'gtexV8.outliers.v8ciseQTLs.globalOutliers.removed.medz.enloc.txt'))
medz_top = filter(enloc_sum,rcp>0.5)
medz_top_enloc_traits = medz_top %>% group_by(trait) %>% 
  mutate(NColocs = length(unique(Gene))) %>%
  sample_n(1) %>% ungroup() %>% filter(NColocs<7) %>%
  mutate(Method='enloc')

coloc_sum = fread(paste0(data_dir, 'gtexV8.outliers.v8ciseQTLs.globalOutliers.removed.medz.coloc.txt'))
medz_top = filter(coloc_sum,p4>0.5)
medz_top_coloc_traits = medz_top %>% group_by(trait) %>% 
  mutate(NColocs = length(unique(Gene))) %>%
  sample_n(1) %>% ungroup() %>% filter(NColocs<7) %>%
  mutate(Method='coloc')

# both_traits = rbind(medz_top_enloc_traits %>% select(trait,NColocs,Method),
#                     medz_top_coloc_traits %>% select(trait,NColocs,Method)) %>% 
#   arrange(by=-NColocs)
# both_traits$trait = factor(both_traits$trait,levels=unique(both_traits$trait))
# ggplot(both_traits, aes(x=trait,y=NColocs, Group=Method)) +
#   ylab('Number genes with coloc prob > 0.5') +
#   geom_bar(stat='identity',aes(fill=Method)) + 
#   theme_bw() + coord_flip() + xlab('') +
#   theme(axis.title = element_text(size=20),
#         axis.text.x = element_text(size=20),
#         axis.text.y = element_text(size=14),
#         title = element_text(size=20),
#         legend.text=element_text(size=20),
#         legend.title=element_text(size=20)) +
#   theme(panel.border = element_blank(),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"))
#   
# gtex_meta = fread('/Users/nicoleferraro/durga_local/data/goats_data/v8_data/GTEx_Analysis_2017-06-05_v8_Annotations_SubjectPhenotypesDS.txt')

# Checking distribution of AFs of relevant variants
gdir = Sys.getenv('RAREDIR')
colnames(enloc_sum)[5] = 'prob'
colnames(coloc_sum)[5] = 'prob'
all_rare_SNPs = fread(paste0(gdir, '/features_v8/byGene/10kb_genebody/all_rare_variants_SNPs_10kb_genebody.txt'),header=F)
colnames(all_rare_SNPs) = c('Ind','Gene','Chr','Start','End')
enloc_with_rare = merge(enloc_sum, all_rare_SNPs, by=c('Ind','Gene'),all.x=T) %>%
  mutate(Method='enloc')
coloc_with_rare = merge(coloc_sum,all_rare_SNPs,by=c('Ind','Gene'),all.x=T) %>%
  mutate(Method='coloc')

coloc_with_rare = rbind(coloc_with_rare,enloc_with_rare)
coloc_with_rare$Chr = sapply(coloc_with_rare$Chr, function(x) strsplit(x,'chr')[[1]][2])
gnomad_afs = fread(paste0(gdir, '/features_v8/gnomad/SNPS.all.INFO'))
colnames(gnomad_afs)[1:2] = c('Chr','Start')
coloc_with_rare = merge(coloc_with_rare,gnomad_afs,by=c('Chr','Start'))
write.table(coloc_with_rare,file=paste0(data_dir, 'gtexV8.coloc.with.rare.txt'),sep='\t',quote=F,row.names=F)

coloc_with_rare$AF_NFE = as.numeric(coloc_with_rare$AF_NFE)
coloc_with_rare = coloc_with_rare %>% mutate(UID=paste(Chr,Start,sep='_')) %>%
  filter(prob>0.5)
coloc_sum = coloc_with_rare %>% group_by(UID) %>% top_n(1,prob) %>% 
  sample_n(1) %>% ungroup()

ggplot(coloc_sum,aes(AF_NFE)) + geom_histogram() + theme_bw() +
    theme(axis.title = element_text(size=20),
          axis.text.x = element_text(size=18),
          axis.text.y = element_text(size=18),
          title = element_text(size=20)) +
    theme(panel.border = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"))

coloc_traits = coloc_with_rare %>% group_by(trait,Gene) %>% sample_n(1) %>% ungroup()
trait_table = as.data.frame(table(coloc_traits$trait))
ggplot(trait_table,aes(Freq)) + geom_histogram(bins=10) + theme_bw() +
  ylab('Number of traits') + xlab('Number of genes') +
  theme(axis.title = element_text(size=20),
        axis.text.x = element_text(size=18),
        axis.text.y = element_text(size=18),
        title = element_text(size=20)) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

### sqtl analysis
data_dir = '/Users/nicoleferraro/durga_local/data/goats_data/v8_data/splicing/'
sqtl_results = fread(paste0(data_dir, 'predixcan_splicing_significant.csv'))
sqtl_results$Start = sapply(sqtl_results$gene, function(x) strsplit(x,'_')[[1]][3])
sqtl_results$End = sapply(sqtl_results$gene, function(x) strsplit(x,'_')[[1]][4])
sqtl_results$Start = as.numeric(sqtl_results$Start)
sqtl_results$End = as.numeric(sqtl_results$End)
gene_positions = fread('/Users/nicoleferraro/durga_local/data/goats_data/v8_data/gencode.v26.GRCh38.genes.gtf')

splicing_variants = fread(paste0(data_dir, 'gtexV8_splicing_outliers_rare_snps_position.txt'))
splicing_variants$chr = sapply(splicing_variants$Chr, function(x) strsplit(x,'chr')[[1]][2])

get_sqtl_overlap <- function(i) {
  schr = splicing_variants$chr[i]
  spos = splicing_variants$Pos[i]
  sqtl_chr = filter(sqtl_results,chr == schr,Start-10000 < spos && End+10000 > spos)
  if (nrow(sqtl_chr) > 0) {
    sqtl_chr$Gene = splicing_variants$Gene[i]
    sqtl_chr$VarChr = schr
    sqtl_chr$VarPos = spos
    return(sqtl_chr %>% select(zscore,pvalue,best_gwas_p,tissue,trait,pos1,pos2,Start,End,Gene,VarChr,VarPos))
  }
}

all_overlap = do.call(rbind, lapply(1:nrow(splicing_variants),get_sqtl_overlap))


