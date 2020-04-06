#!/usr/bin/env Rscript

## Plot enrichment of single tissue z-scores

library(RColorBrewer)
library(data.table)
library(dplyr)
library(ggplot2)
library(epitools)

data_dir = '/users/nferraro/data/goats_data/v8_data/'

medz_single = fread(paste0(data_dir, 'tissue_specific_outliers/gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.medz.top.single.filt20.txt'),data.table=F)
medz_outliers_only = filter(medz_single,abs(value)>2) %>% mutate(Y='outlier') %>%
  mutate(UID=paste(variable,Gene,sep='_'))
medz_controls_only = filter(medz_single,abs(value)<2,Gene %in% medz_outliers_only$Gene) %>% 
  mutate(UID=paste(variable,Gene,sep='_')) %>% mutate(Y='control')

rare_SVs = fread(paste0(data_dir, 'all_genes_with_rare_HallLabSV.txt'),header=F) %>%
  mutate(UID = paste(V1,V2,sep='_'))
rare_SNPs = fread(paste0(data_dir, 'all_genes_with_rare_SNPs.txt'),header=F) %>%
  mutate(UID = paste(V1,V2,sep='_')) %>% filter(V1 %in% rare_SVs$V1)
rare_indels = fread(paste0(data_dir, 'all_genes_with_rare_indels.txt'),header=F) %>%
  mutate(UID = paste(V1,V2,sep='_')) %>% filter(V1 %in% rare_SVs$V1)

medz_outliers_only = filter(medz_outliers_only, variable %in% rare_SVs$V1)
medz_controls_only = filter(medz_controls_only, variable %in% rare_SVs$V1)

zts = c(3,4,5,6,7,8,9,10)
counts = data.frame(Over = numeric(), Under = numeric(), ZT = numeric())
all_tissue_table = data.frame(Tissue = character(), NumOuts = numeric(), ZT = numeric())
for (zt in zts) {
  no = nrow(filter(medz_outliers_only,value > zt))
  nu = nrow(filter(medz_outliers_only,value < -zt))
  counts = rbind(counts, data.frame(Over = no, Under = nu, ZT = zt))
  tissue_table = as.data.frame(table(filter(medz_outliers_only,abs(value)>zt)$Tissue))
  colnames(tissue_table) = c('Tissue','NumOuts')
  tissue_table$ZT = zt
  all_tissue_table = rbind(all_tissue_table, tissue_table)
}

counts.melted = melt(counts, id.vars='ZT')
counts = counts %>% mutate(Ratio = Over/Under)
ggplot(counts,aes(x=ZT,y=Ratio)) +
  geom_bar(stat='identity') +
  scale_fill_manual(values=c('darkgrey','lightgrey')) +
  theme_bw() + ylab('Ratio of over:under expression') + xlab('Z-score threshold') +
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=20),
        title = element_text(size=20),
        legend.text=element_text(size=20),
        legend.title=element_blank()) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggplot(counts.melted,aes(x=ZT,y=value,Group=variable)) +
  geom_bar(stat='identity',aes(fill=variable),position='dodge') +
  scale_fill_manual(values=c('darkgrey','lightgrey')) +
  theme_bw() + ylab('Number of outliers') + xlab('Z-score threshold') +
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=20),
        title = element_text(size=20),
        legend.text=element_text(size=20),
        legend.title=element_blank()) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

samp_nums = fread(paste0(data_dir, 'sample_info/gtexV8_tissue_sample_numbers.txt')) %>%
  arrange(by=Freq)
all_tissue_table$ZT = factor(all_tissue_table$ZT, levels=unique(all_tissue_table$ZT))
all_tissue_table$Tissue = factor(all_tissue_table$Tissue, levels=samp_nums$Var1)
ggplot(all_tissue_table,aes(x=Tissue,y=NumOuts,Group=ZT)) +
  geom_bar(stat='identity',aes(fill=ZT)) +
  theme_bw() + ylab('Number of outliers') + xlab('Tissue') +
  theme(axis.title = element_text(size=20),
        axis.text = element_blank(),
        title = element_text(size=20),
        legend.text=element_text(size=20),
        legend.title=element_blank()) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

risks = data.frame(riskratio = numeric(), lower = numeric(), upper = numeric(), p.value = numeric(), Method = character(), Type = character(), ZT = numeric(), NOutliers = numeric())
for (zt in zts) {
  print(zt)
  moutliers = filter(medz_outliers_only) %>% filter(abs(value) > zt)
  # moutliers = filter(medz_single) %>% group_by(Gene,Tissue) %>%
  #   top_n(1,abs(value)) %>% ungroup() %>% filter(value > zt)
  mcount = nrow(moutliers)
  mcontrols = filter(medz_controls_only, Gene %in% moutliers$Gene)
  #mcontrols = filter(medz_single, Gene %in% moutliers$Gene, !(UID %in% moutliers$UID))

  z_snp_table = rbind(c(length(which(!(mcontrols$UID %in% rare_SNPs$UID))),
                        length(which(mcontrols$UID %in% rare_SNPs$UID))),
                      c(length(which(!(moutliers$UID %in% rare_SNPs$UID))),
                        length(which(moutliers$UID %in% rare_SNPs$UID))))
  z_snp_rr = epitab(z_snp_table,method="riskratio")$tab
  z_indel_table = rbind(c(length(which(!(mcontrols$UID %in% rare_indels$UID))),
                          length(which(mcontrols$UID %in% rare_indels$UID))),
                        c(length(which(!(moutliers$UID %in% rare_indels$UID))),
                          length(which(moutliers$UID %in% rare_indels$UID))))
  z_indel_rr = epitab(z_indel_table,method="riskratio")$tab
  z_sv_table = rbind(c(length(which(!(mcontrols$UID %in% rare_SVs$UID))),
                       length(which(mcontrols$UID %in% rare_SVs$UID))),
                     c(length(which(!(moutliers$UID %in% rare_SVs$UID))),
                       length(which(moutliers$UID %in% rare_SVs$UID))))
  z_sv_rr = epitab(z_sv_table,method="riskratio")$tab
  
  rr_table = rbind(c(z_snp_rr[2,5], z_snp_rr[2,6], z_snp_rr[2,7], z_snp_rr[2,8], 'Zscore', 'SNPs', zt, mcount),
                   c(z_indel_rr[2,5], z_indel_rr[2,6], z_indel_rr[2,7], z_indel_rr[2,8], 'Zscore', 'indels', zt, mcount),
                   c(z_sv_rr[2,5], z_sv_rr[2,6], z_sv_rr[2,7], z_sv_rr[2,8], 'Zscore', 'SVs', zt, mcount))
  risks = rbind(risks, rr_table)
}

colnames(risks) = c('Riskratio', 'Lower', 'Upper', 'Pval', 'Method', 'Type', 'ZT', 'NOutliers')
risks$Riskratio = as.numeric(as.character(risks$Riskratio))
risks$Lower = as.numeric(as.character(risks$Lower))
risks$Upper = as.numeric(as.character(risks$Upper))

risks_over = risks %>% mutate(Category = 'Over')
risks = rbind(risks %>% mutate(Category = 'Under'), risks_over)

ggplot(filter(risks, Type == 'indels'), aes(x=ZT,y=Riskratio)) +
  geom_point(size=4) + theme_bw() +
  geom_hline(yintercept=1,color='black',linetype='dashed') +
  scale_color_manual(values=c('darkgrey','lightgrey')) +
  ggtitle('SVs') +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), position=position_dodge(), width=0) +
  xlab('Z-score threshold') +
  ylab('Relative risk') +
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=20),
        title = element_text(size=20),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20),
        legend.key.size = unit(1.5, 'lines')) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

