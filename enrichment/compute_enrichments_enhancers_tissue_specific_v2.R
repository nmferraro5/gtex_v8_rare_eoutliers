library(data.table)
library(dplyr)
library(ggplot2)
library(epitools)

data_dir = '/users/nferraro/data/goats_data/v8_data/'

all_enh = fread(paste0(data_dir, 'enhancers/Enh_links_all.txt'),data.table=F)
all_enh = all_enh %>% mutate(UID = paste(Chr,Start,End,Gene,sep='_')) %>%
  group_by(UID) %>% mutate(NTissues = length(unique(Tissue))) %>% ungroup()
tissue.map = fread(paste0(data_dir, 'roadmap_gtex_map_castel.txt'),data.table=F) %>%
  select(GTEx_tissue,EID,Tissue.Name) 
colnames(all_enh)[7] = 'EID'
all_enh = merge(all_enh, tissue.map, by='EID', all.x=T) %>%
  filter(!is.na(GTEx_tissue))
keep_enh = fread(paste0(data_dir, 'enhancers/enh_keep_n1.txt'))
colnames(keep_enh) = c('Chr', 'Start', 'End', 'Gene', 'EID')
all_enh = all_enh %>% mutate(UID = paste(Gene,Ind,GTEx_tissue,sep='_'))
keep_enh = merge(keep_enh, all_enh, by=c('Chr','Start','End','Gene', 'EID'))
rm(all_enh)
gtex.map = fread(paste0(data_dir, 'gtex.tissue.map.txt'))
colnames(gtex.map)[2] = 'GTEx_tissue'
keep_enh = merge(keep_enh, gtex.map %>% select(BroadName,GTEx_tissue), by='GTEx_tissue')
colnames(keep_enh)[13] = 'BroadTissue'
keep_enh = keep_enh %>% mutate(UID = paste(Gene,Ind,BroadTissue,sep='_'))
keep_enh = filter(keep_enh, VarType != 'HallLabSV', AF < 0.005) # worked well with FDR < 0.0027 and no SVS and AF < 0.005
keep_enh = keep_enh %>% mutate(GID = paste(Gene,Ind,sep='_'))

ts.outliers = fread(paste0(data_dir, 'tissue_specific_outliers/gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.subset.knn.FDR.specific.filtered.txt'))
ts.outliers = ts.outliers %>% mutate(UID = paste(Gene,Ind,sep='_'))
ts.outliers = ts.outliers %>% group_by(UID) %>% top_n(1,abs(value)) %>% ungroup()
ts.controls = fread(paste0(data_dir, 'tissue_specific_outliers/gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.subset.knn.FDR.specific.top5.controls.txt'))
ts.controls$FDR = rep(NA,nrow(ts.controls))
ts.outliers = filter(ts.outliers, BroadTissue %in% keep_enh$BroadTissue)
ts.outliers = filter(ts.outliers, FDR < 0.0027)
ts.outliers$Gene = sapply(ts.outliers$Gene, function(x) strsplit(x, '[.]')[[1]][1])
ts.controls$Gene = sapply(ts.controls$Gene, function(x) strsplit(x, '[.]')[[1]][1])
ts.outliers = filter(ts.outliers, Gene %in% keep_enh$Gene, Ind %in% keep_enh$Ind)
ts.controls = filter(ts.controls, Gene %in% keep_enh$Gene, Gene %in% ts.outliers$Gene, Ind %in% keep_enh$Ind)
ts.controls = merge(ts.controls, ts.outliers %>% select(Gene,BroadTissue), by='Gene')
ts.outliers = ts.outliers %>% mutate(UID = paste(Gene,Ind,BroadTissue,sep='_')) %>% mutate(GID = paste(Gene,Ind,sep='_'))
ts.controls = ts.controls %>% mutate(UID = paste(Gene,Ind,BroadTissue,sep='_')) %>% mutate(GID = paste(Gene,Ind,sep='_'))

print('Starting relative risk calc')
risks = data.frame(riskratio = numeric(), lower = numeric(), upper = numeric(), p.value = numeric(), ZT = numeric(), NOutliers = numeric(), Type = character())
fdr_thresh = c(0.05, 0.01, 0.001, 0.0001, 0.00001, 0.000001)
zt_thresh = c(3,5,7,9)
for (zt in zt_thresh) {
  print(zt)
  #filt.outliers = ts.outliers %>% filter(FDR < ft)
  filt.outliers = ts.outliers %>% filter(abs(value)>zt)
  mcount = nrow(filt.outliers)
  filt.controls = filter(ts.controls, Gene %in% filt.outliers$Gene, BroadTissue %in% filt.outliers$BroadTissue)
  no_no = nrow(filter(filt.controls, !(UID %in% keep_enh$UID)))
  no_yes = nrow(filter(filt.controls, UID %in% keep_enh$UID))
  yes_no = nrow(filter(filt.outliers, !(UID %in% keep_enh$UID)))
  yes_yes = nrow(filter(filt.outliers, UID %in% keep_enh$UID))
  
  counttable = rbind(c(no_no,no_yes),c(yes_no,yes_yes))
  print(counttable)
  enh_rr = epitab(counttable,method="riskratio")$tab
  risks = rbind(risks,c(enh_rr[2,5], enh_rr[2,6], enh_rr[2,7], enh_rr[2,8], zt, mcount, 'Matched'))
}

colnames(risks) = c('Riskratio', 'Lower', 'Upper', 'Pval', 'ZT', 'Count', 'Type')
#risks$FT = factor(risks$FT, levels=fdr_thresh)
write.table(risks, file=paste0(data_dir, 'tissue_specific_outliers/gtexV8.cor.FDR.tissue.specific.enhancer.risks.txt'),sep='\t',quote=F,row.names=F)

# risks_all$AF = factor(risks_all$AF, levels=c(0.01, 0.001))
# ggplot(filter(risks_all, Type != 'SVs'), aes(x=FT,y=Riskratio,Group=AF)) +
#   geom_point(size=7,aes(color=AF),position=position_dodge(width=0.2)) + theme_bw() +
#   geom_errorbar(aes(ymin = Lower, ymax = Upper), position=position_dodge(width=0.2), width=0) +
#   geom_hline(yintercept=1) + xlab('FDR threshold') +
#   scale_color_manual(values=c('darkgrey', 'black')) +
#   theme(axis.title = element_text(size=24),
#         axis.text = element_text(size=22),
#         title = element_text(size=20),
#         legend.text=element_text(size=22),
#         legend.title=element_text(size=22),
#         legend.key.size = unit(1.5, 'lines')) +
#   theme(panel.border = element_blank(),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"))

