## Assess tissue-specific enhancer enrichments

library(data.table)
library(dplyr)
library(ggplot2)
require(foreach)
require(doMC)
library(gplots)
require(RColorBrewer)
library(epitools)

data_dir = '/users/nferraro/data/goats_data/v8_data/'
#outliers.controls.specific = fread(paste0(data_dir, 'tissue_specific_outliers/gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.subset.knn.FDR.specific.txt'))
#outliers = filter(outliers.controls.specific, !is.na(Tissue)) %>%
#  filter(FDR < 0.0027) %>% mutate(Y=1)
#controls = filter(outliers.controls.specific, is.na(Tissue), Gene %in% outliers$Gene) %>%
#  filter(is.na(value) | abs(value) < 3) %>% mutate(Y=0)
#controls$Tissue = sample(outliers$Tissue, nrow(controls), replace=T)
#outliers.controls.specific = rbind(outliers,controls)
outliers.controls.specific = fread(paste0(data_dir, 'gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.medz.filtered.txt'))
outliers = filter(outliers.controls.specific, abs(MedZ) > 3)
controls = filter(outliers.controls.specific, abs(MedZ) < 3, Gene %in% outliers$Gene)
outliers.controls.specific = rbind(outliers,controls)


all_enh = fread(paste0(data_dir, 'enhancers/Enh_links_all.txt'),data.table=F)
tissue.map = fread(paste0(data_dir, 'roadmap_gtex_map_castel.txt'),data.table=F) %>%
  select(GTEx_tissue,EID,Tissue.Name) 
colnames(all_enh)[7] = 'EID'
all_enh = merge(all_enh, tissue.map, by='EID', all.x=T) %>%
  filter(!is.na(GTEx_tissue))

keep_enh = fread(paste0(data_dir, 'enhancers/enh_keep_n1.txt'))
colnames(keep_enh) = c('Chr', 'Start', 'End', 'Gene', 'EID')
keep_enh = merge(keep_enh, all_enh, by=c('Chr','Start','End','Gene', 'EID'))
rm(all_enh)

kids = unique(filter(keep_enh, VarType == 'HallLabSV')$Ind)

outliers.controls.specific$GeneID = sapply(outliers.controls.specific$Gene, function(x) strsplit(x,'[.]')[[1]][1])
outliers.controls.specific = outliers.controls.specific %>% filter(Ind %in% kids)
colnames(keep_enh)[4] = 'GeneID'
keep_enh = keep_enh %>% select(GeneID,Ind,AF,VarType,GTEx_tissue) %>% mutate(UID = paste0(GeneID,Ind))
outliers.controls.specific$GeneID = sapply(outliers.controls.specific$Gene, function(x) strsplit(x, '[.]')[[1]][1])
outliers.controls.enh = merge(outliers.controls.specific, keep_enh, by=c('GeneID','Ind'), all.x=T)

#outliers.controls.enh = outliers.controls.enh %>% group_by(GeneID,Ind,Tissue) %>% 
#  mutate(HasRV = ifelse(any(GTEx_tissue == Tissue), 1, 0)) %>% 
#  mutate(NT = length(unique(GTEx_tissue))) %>% ungroup()
#
#risks = data.frame(Risk = numeric(), Lower = numeric(), Upper = numeric(), Pval = numeric(), ZT = numeric(), NOutliers = numeric())
#zts = c(3,4,5,6)
#for (zt in zts) {
#  print(paste0('ZT = ', zt))
#  outliers = filter(outliers.controls.enh, Y == 1, abs(value) > zt)
#  controls = filter(outliers.controls.enh, Y == 0, Gene %in% outliers$Gene)
#  no_no = length(unique(filter(controls, HasRV == 0)$UID))
#  no_yes = length(unique(filter(controls, HasRV == 1)$UID))
#  yes_no = length(unique(filter(outliers, HasRV == 0)$UID))
#  yes_yes = length(unique(filter(outliers, HasRV == 1)$UID))
#  counttable = rbind(c(no_no, no_yes), c(yes_no, yes_yes))
#  rr = epitab(counttable, method = 'riskratio')
#  risks = rbind(risks, data.frame(Risk = rr$tab[2,5],
#                                  Lower = rr$tab[2,6],
#                                  Upper = rr$tab[2,7],
#                                  Pval = rr$tab[2,8],
#                                  ZT = zt,
#                                  NOutliers = mcount))
#}
###

#outliers.controls.enh = merge(outliers.controls.specific, si_enh, by=c('GeneID','Ind'),all.x=T)
#outliers.controls.enh = outliers.controls.enh %>% select(-Gene,-NTissues,-Tissue.Name,-EID)

#outliers.controls.enh = outliers.controls.enh %>% group_by(GeneID,Ind) %>%
#  mutate(NTissues = length(unique(which(!is.na(GTEx_tissue))))) %>%
#  ungroup() %>% mutate(Match = ifelse(any(Tissue == GTEx_tissue), 1, 0)) %>%
#  mutate(UID = paste0(GeneID,Ind)) 
#tkeep = unique(outliers.controls.enh$GTEx_tissue)
#outliers.controls.enh = filter(outliers.controls.enh, Tissue %in% tkeep)
outliers.controls.enh = outliers.controls.enh %>% group_by(GeneID,Ind) %>%
  mutate(NTissues = length(unique(which(!is.na(GTEx_tissue))))) %>%
  ungroup() %>% mutate(UID = paste0(GeneID,Ind)) 
outliers.controls.specific = outliers.controls.specific %>% mutate(UID = paste0(GeneID,Ind))

risks = data.frame(Risk = numeric(), Lower = numeric(), Upper = numeric(), Pval = numeric(), ZT = numeric(), NOutliers = numeric())
zts = c(3,4,5,6)
for (zt in zts) {
  print(paste0('ZT = ', zt))
  ft = 2*pnorm(-abs(zt))
  if (zt > 0) {
    #outliers.specific = filter(outliers.controls.specific,!is.na(Tissue)) %>%
    #  filter(FDR < ft, abs(value) > 3)
    outliers.specific = filter(outliers.controls.specific,abs(MedZ) > zt)
  } else {
    outliers.specific = filter(outliers.controls.specific,!is.na(Tissue)) %>%
      filter(value < zt) 
  }
  mcount = nrow(outliers.specific)
  outliers.controls = filter(outliers.controls.specific, abs(MedZ) < 3) %>%
    filter(Gene %in% outliers.specific$Gene)
  
  no_no = nrow(filter(outliers.controls, !(UID %in% keep_enh$UID)))
  no_yes = nrow(filter(outliers.controls, UID %in% keep_enh$UID))
  yes_no = nrow(filter(outliers.specific, !(UID %in% keep_enh$UID)))
  yes_yes = nrow(filter(outliers.specific, UID %in% keep_enh$UID))
  counttable = rbind(c(no_no, no_yes), c(yes_no, yes_yes))
  rr = epitab(counttable, method = 'riskratio')
  risks = rbind(risks, data.frame(Risk = rr$tab[2,5],
                                  Lower = rr$tab[2,6],
                                  Upper = rr$tab[2,7],
                                  Pval = rr$tab[2,8],
                                  ZT = zt,
                                  NOutliers = mcount))
}

write.table(risks, file=paste0(data_dir, 'tissue_specific_outliers/medz_gene_enhancer_risks.txt'),sep='\t',quote=F,row.names=F)
