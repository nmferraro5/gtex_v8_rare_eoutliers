#!/usr/bin/env Rscript

# Single tissue outliers

## Get path to the GOATs data directories
dir = Sys.getenv('RAREDIR')

## Load required packages
require(data.table)
require(dplyr)
library(ggplot2)
library(epitools)
library(RColorBrewer)

data_dir = '/Users/nicoleferraro/durga_local/data/goats_data/v8_data/'
singleZ_data = fread(paste0(data_dir, 'gtexV8.outliers.controls.single.tissue.exp.txt'))
singleZ_outliers = filter(singleZ_data, abs(value) > 3)

cor_splice = fread('/Users/nicoleferraro/durga_local/data/goats_data/v8_data/splicing/gtexV8.splice.outlier.controls.v8ciseQTLs.globalOutliers.removed.subset.knn.txt')
cor_outliers = filter(cor_splice,FDR < 1e-09) %>% mutate(Y='outlier')
cor_controls = filter(cor_splice,Gene %in% cor_outliers$Gene,FDR > 1e-09) %>% mutate(Y='control')
cor_splice = rbind(cor_outliers,cor_controls) %>% mutate(UID=paste(Ind,Gene,sep='_'))

single_splice = fread(paste0(data_dir, 'splicing/gtexV8.outliers.controls.tissue.specific.splicing.txt'))
pthresh = 2*pnorm(-abs(3))
single_outliers = filter(single_splice,!is.na(TissueName),value < pthresh) %>%
  mutate(Y='outlier') %>% top_n(nrow(cor_outliers),-value)
single_controls = filter(single_splice,is.na(TissueName),Gene %in% single_outliers$Gene) %>%
  mutate(Y='control')
single_splice = rbind(single_outliers,single_controls) %>%
  mutate(UID=paste(Ind,Gene,sep='_'))

rare_snps = fread(paste0(data_dir, 'all_genes_with_rare_SNPs.txt'),header=F) %>%
  mutate(UID=paste(V1,V2,sep='_'))
rare_indels = fread(paste0(data_dir, 'all_genes_with_rare_indels.txt'),header=F) %>%
  mutate(UID=paste(V1,V2,sep='_'))
rare_SVs = fread(paste0(data_dir, 'all_genes_with_rare_HallLabSV.txt'),header=F) %>%
  mutate(UID=paste(V1,V2,sep='_'))

single_splice = single_splice %>% mutate(HasSNP = ifelse(UID %in% rare_snps$UID, 1, 0)) %>%
  mutate(HasIndel = ifelse(UID %in% rare_indels$UID, 1, 0)) %>%
  mutate(HasSV = ifelse(UID %in% rare_SVs$UID, 1, 0))

cor_splice = cor_splice %>% mutate(HasSNP = ifelse(UID %in% rare_snps$UID, 1, 0)) %>%
  mutate(HasIndel = ifelse(UID %in% rare_indels$UID, 1, 0)) %>%
  mutate(HasSV = ifelse(UID %in% rare_SVs$UID, 1, 0))

snp_table = rbind(table(filter(single_splice,Y=='control')$HasSNP),table(filter(single_splice,Y=='outlier')$HasSNP))
snp_rr = epitab(snp_table,method="riskratio")$tab
print(snp_rr)

indel_table = rbind(table(filter(single_splice,Y=='control')$HasIndel),table(filter(single_splice,Y=='outlier')$HasIndel))
indel_rr = epitab(indel_table,method="riskratio")$tab
print(indel_rr)

sv_table = rbind(table(filter(single_splice,Y=='control')$HasSV),table(filter(single_splice,Y=='outlier')$HasSV))
sv_rr = epitab(sv_table,method="riskratio")$tab
print(sv_rr)

all_outliers_rr = rbind(c(snp_rr[2,4],snp_rr[2,5],snp_rr[2,6],snp_rr[2,7],'SNPs'),
               c(indel_rr[2,4],indel_rr[2,5],indel_rr[2,6],indel_rr[2,7],'Indels'),
               c(sv_rr[2,4],sv_rr[2,5],sv_rr[2,6],sv_rr[2,7],'SVs'))

top_outliers_rr = rbind(c(snp_rr[2,5],snp_rr[2,6],snp_rr[2,7],snp_rr[2,8],'SNPs'),
                        c(indel_rr[2,5],indel_rr[2,6],indel_rr[2,7],indel_rr[2,8],'Indels'),
                        c(sv_rr[2,5],sv_rr[2,6],sv_rr[2,7],sv_rr[2,8],'SVs'))

snp_table = rbind(table(filter(cor_splice,Y=='control')$HasSNP),table(filter(cor_splice,Y=='outlier')$HasSNP))
snp_rr = epitab(snp_table,method="riskratio")$tab
print(snp_rr)

indel_table = rbind(table(filter(cor_splice,Y=='control')$HasIndel),table(filter(cor_splice,Y=='outlier')$HasIndel))
indel_rr = epitab(indel_table,method="riskratio")$tab
print(indel_rr)

sv_table = rbind(table(filter(cor_splice,Y=='control')$HasSV),table(filter(cor_splice,Y=='outlier')$HasSV))
sv_rr = epitab(sv_table,method="riskratio")$tab
print(sv_rr)

cor_outliers_rr = rbind(c(snp_rr[2,5],snp_rr[2,6],snp_rr[2,7],snp_rr[2,8],'SNPs'),
                        c(indel_rr[2,5],indel_rr[2,6],indel_rr[2,7],indel_rr[2,8],'Indels'),
                        c(sv_rr[2,5],sv_rr[2,6],sv_rr[2,7],sv_rr[2,8],'SVs'))

all_rr = rbind(as.data.frame(top_outliers_rr) %>% mutate(Method = 'Single'),
               as.data.frame(cor_outliers_rr) %>% mutate(Method = 'COR'))

all_rr$riskratio = as.numeric(as.character(all_rr$riskratio))
all_rr$lower = as.numeric(as.character(all_rr$lower))
all_rr$upper = as.numeric(as.character(all_rr$upper))

meth.colors = brewer.pal(9, 'YlGnBu')[3:9]
plotcols = c(meth.colors[2],'grey')
names(plotcols) = c('COR','Single')
ggplot(all_rr, aes(x=V5,y=riskratio,Group=Method)) +
  geom_point(size = 6, position=position_dodge(width = 0.7),aes(color=Method)) +
  geom_errorbar(aes(ymin = lower, ymax = upper), position=position_dodge(width = 0.7), width=0) +
  theme_bw() + ylab('Relative risk') + xlab('') +
  ggtitle('') + scale_color_manual(values=plotcols) +
  geom_hline(yintercept=1, linetype="dashed", colour="darkgrey") +
  theme(axis.title = element_text(size=20),
        axis.text.y = element_text(size=20),
        axis.text.x = element_text(size=20,hjust=1,angle=25),
        title = element_text(size=20),
        legend.text=element_text(size=18)) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))






