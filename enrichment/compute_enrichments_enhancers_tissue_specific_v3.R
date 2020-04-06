library(data.table)
library(dplyr)
library(epitools)
library(ggplot2)

data_dir = '/users/nferraro/data/goats_data/v8_data/'

cor_ts = fread(paste0(data_dir, 'tissue_specific_outliers/gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.subset.knn.FDR.specific.top5.txt'))

all_enh = fread(paste0(data_dir, 'enhancers/Enh_links_all.txt'),data.table=F)
tissue.map = fread(paste0(data_dir, 'roadmap_gtex_map_castel.txt'),data.table=F) %>%
  select(GTEx_tissue,EID,Tissue.Name) 
colnames(all_enh)[7] = 'EID'
all_enh = merge(all_enh, tissue.map, by='EID', all.x=T) %>%
  filter(!is.na(GTEx_tissue))

keep_enh = fread(paste0(data_dir, 'enhancers/enh_keep_n1.txt'))
colnames(keep_enh) = c('Chr', 'Start', 'End', 'Gene', 'EID')
keep_enh = merge(keep_enh, all_enh, by=c('Chr','Start','End','Gene', 'EID'))
colnames(keep_enh)[9] = 'Tissue'
keep_enh = filter(keep_enh, VarType != 'HallLabSV')
rm(all_enh)

cor_ts$Gene = sapply(cor_ts$Gene, function(x) strsplit(x, '[.]')[[1]][1])
cor_ts_outliers = filter(cor_ts, !is.na(FDR)) 
cor_ts_single = filter(cor_ts_outliers, NT == 1) %>%
  group_by(Ind,Gene) %>% top_n(2,abs(value)) %>%
  mutate(MinVal = min(abs(value))) %>% mutate(MaxVal = max(abs(value))) %>%
  mutate(IsSingle = ifelse(MinVal < 2 & MaxVal > 3, 1, 0)) %>%
  filter(IsSingle == 1) %>% top_n(1,abs(value)) %>% ungroup() %>% 
  filter(Tissue %in% keep_enh$Tissue, Ind %in% keep_enh$Ind) %>%
  select(-IsSingle,-MinVal,-MaxVal)

risks = data.frame(Risk = numeric(), Lower = numeric(), Upper = numeric(), Pval = numeric(), ZT = numeric(), Type = character())
zts = c(3,4,5,6)
for (zt in zts) {
  ft = 2*pnorm(-abs(zt))
  ts_enh_data = merge(filter(cor_ts_single, FDR < ft), keep_enh %>% select(-Tissue.Name), by=c('Ind','Gene'),all.x=T) %>%
    mutate(HasMatch = ifelse(Tissue.x == Tissue.y, 1, 0)) %>%
    mutate(HasMatch = ifelse(is.na(HasMatch), 0, HasMatch)) %>%
    mutate(HasOther = ifelse(!is.na(Tissue.y) & Tissue.x != Tissue.y, 1, 0)) %>%
    mutate(HasOther = ifelse(is.na(HasOther), 0, HasOther)) %>%
    group_by(Gene,Ind,Tissue.x) %>% mutate(Match = max(HasMatch,na.rm=T)) %>%
    mutate(Other = max(HasOther, na.rm=T)) %>% sample_n(1) %>% ungroup() %>%
    select(-HasMatch, -HasOther)
  ts_controls = filter(cor_ts, is.na(FDR), Gene %in% ts_enh_data$Gene)
  #ts_controls$Tissue = sample(cor_ts$Tissue, nrow(ts_controls), replace=T)
  ts_controls$Tissue = sapply(ts_controls$Gene, function(x) ts_enh_data$Tissue.x[sample(which(ts_enh_data$Gene == x),1)])
  ts_controls = merge(ts_controls, keep_enh %>% select(-Tissue.Name), by=c('Ind','Gene'),all.x=T) %>%
    mutate(HasMatch = ifelse(Tissue.x == Tissue.y, 1, 0)) %>%
    mutate(HasMatch = ifelse(is.na(HasMatch), 0, HasMatch)) %>%
    mutate(HasOther = ifelse(!is.na(Tissue.y) & Tissue.x != Tissue.y, 1, 0)) %>%
    mutate(HasOther = ifelse(is.na(HasOther), 0, HasOther)) %>%
    group_by(Gene,Ind,Tissue.x) %>% mutate(Match = max(HasMatch,na.rm=T)) %>%
    mutate(Other = max(HasOther, na.rm=T)) %>% sample_n(1) %>% ungroup() %>%
    select(-HasMatch, -HasOther)
  
  no_no = nrow(filter(ts_controls, Match == 0))
  no_yes = nrow(filter(ts_controls, Match == 1))
  yes_no = nrow(filter(ts_enh_data, Match == 0))
  yes_yes = nrow(filter(ts_enh_data, Match == 1))
  matchtable = rbind(c(no_no,no_yes),c(yes_no,yes_yes))
  rr = epitab(matchtable, method = 'riskratio')
  risks = rbind(risks, data.frame(Risk = rr$tab[2,5],
                                  Lower = rr$tab[2,6],
                                  Upper = rr$tab[2,7],
                                  Pval = rr$tab[2,8],
                                  ZT = zt,
                                  Type = 'Match'))
  
  no_no = nrow(filter(ts_controls, Other == 0))
  no_yes = nrow(filter(ts_controls, Other == 1, Match == 0))
  yes_no = nrow(filter(ts_enh_data, Other == 0))
  yes_yes = nrow(filter(ts_enh_data, Other == 1, Match == 0))
  othertable = rbind(c(no_no,no_yes),c(yes_no,yes_yes))
  rr = epitab(othertable, method = 'riskratio')
  risks = rbind(risks, data.frame(Risk = rr$tab[2,5],
                                  Lower = rr$tab[2,6],
                                  Upper = rr$tab[2,7],
                                  Pval = rr$tab[2,8],
                                  ZT = zt,
                                  Type = 'Other'))
  
}


write.table(risks, file=paste0(data_dir, 'tissue_specific_outliers/gtexV8.cor.FDR.tissue.specific.enhancer.noSV.geneControls.risks.txt'),sep='\t',quote=F,row.names=F)

risks = fread(paste0(data_dir, 'tissue_specific_outliers/gtexV8.cor.FDR.tissue.specific.enhancer.noSV.geneControls.risks.txt'))
risks$ZT = factor(risks$ZT, levels=c(3,4,5,6))
ggplot(risks, aes(x=Type,y=Risk,Group=ZT)) +
  geom_point(size = 6, position=position_dodge(width = 0.7),aes(color=ZT)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), position=position_dodge(width = 0.7), width=0) +
  theme_bw() + geom_hline(yintercept=1, linetype="dashed", colour="darkgrey") +
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=20),
        legend.title = element_text(size=20),
        legend.text = element_text(size=18)) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

### scratch
ts_enh_data = merge(filter(cor_ts_single, FDR < ft), keep_enh %>% select(-Tissue.Name), by=c('Ind','Gene'),all.x=T) %>%
  group_by(Gene,Ind,Tissue.x) %>%
  mutate(NMatch = length(which(Tissue.x == Tissue.y))) %>%
  mutate(HasOther = ifelse(!is.na(Tissue.y) & Tissue.x != Tissue.y, 1, 0)) %>%
  mutate(HasOther = ifelse(is.na(HasOther), 0, HasOther)) %>%
  mutate(Other = max(HasOther, na.rm=T)) %>% sample_n(1) %>% ungroup() %>%
  select(-HasOther)
