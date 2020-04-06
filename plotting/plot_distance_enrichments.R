library(ggplot2)
library(data.table)
library(dplyr)
library(RColorBrewer)

gtex_v8_figure_theme <- function() {
  return(theme(plot.title = element_text(face="plain",size=8), text = element_text(size=8),
               axis.text=element_text(size=7), panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(),panel.background = element_blank(), 
               panel.border = element_blank(),
               axis.line = element_line(colour = "black"), 
               legend.text = element_text(size=7), legend.title = element_text(size=8)))
}

data_dir = '/Users/nicoleferraro/durga_local/data/goats_data/v8_data/'

te_dist_files = dir(paste0(data_dir, 'enrichments/'), '*genebody_Z3_FDR0.01.allTissues.v8ciseQTLs.globalOutliers.removed.medz.RData')
te_dist_files = te_dist_files[which(grepl('enrichmentsDIST_downstream_relativeRisk_v8_vg_', te_dist_files))]

ase_dist_files = dir(paste0(data_dir, 'enrichments/'), '*genebody_Z3_FDR0.01.allTissues.v8ciseQTLs.globalOutliers.removed.ase.RData')
ase_dist_files = ase_dist_files[which(grepl('enrichmentsDIST_downstream_relativeRisk_v8_vg_', ase_dist_files))]

splice_dist_files = dir(paste0(data_dir, 'enrichments/'), '*genebody_Z3_FDR0.01.allTissues.v8ciseQTLs.globalOutliers.removed.splicing.RData')
splice_dist_files = splice_dist_files[which(grepl('enrichmentsDIST_downstream_relativeRisk_v8_vg_', splice_dist_files))]

#splice_dist_files = dir(paste0(data_dir, 'enrichments/'), '*genebody_Z3_FDR0.01.allTissues.v8ciseQTLs.globalOutliers.removed.splicing.RData')
#splice_dist_files = splice_dist_files[which(grepl('enrichmentsDIST_upstream_pseudocount_relativeRisk_', splice_dist_files))]

# te_dist_files = dir(paste0(data_dir, 'enrichments/'), '*genebody_Z3_FDR0.01.allTissues.v8ciseQTLs.globalOutliers.removed.medz.RData')
# te_dist_files = te_dist_files[which(grepl('enrichmentsDIST_upstream', te_dist_files))]
# 
# ase_dist_files = dir(paste0(data_dir, 'enrichments/'), '*genebody_Z3_FDR0.01.allTissues.v8ciseQTLs.globalOutliers.removed.ase.RData')
# ase_dist_files = ase_dist_files[which(grepl('enrichmentsDIST_upstream', ase_dist_files))]
# 
# splice_dist_files = dir(paste0(data_dir, 'enrichments/'), '*genebody_Z3_FDR0.01.allTissues.v8ciseQTLs.globalOutliers.removed.splicing.RData')
# splice_dist_files = splice_dist_files[which(grepl('enrichmentsDIST_upstream', splice_dist_files))]

read_dist_data <- function(dfile) {
  print(dfile)
  if (grepl('10kb_genebody', dfile)) {
    window = 0
  } else if (grepl('1bp_200kb', dfile)) {
    window = 200
  } else if (grepl('200kb_400',dfile)) {
    window = 400
  } else if (grepl('400kb_600',dfile)) {
    window = 600
  } else if (grepl('600kb_800',dfile)) {
    window = 800
  } else if (grepl('800kb_1MB',dfile)) {
    window = 1000
  }
  load(paste0(data_dir, 'enrichments/', dfile))
  all_logits_top$Window = window
  if (grepl('ase',dfile)) {
    type = 'ASE'
  } else if (grepl('splicing', dfile)) {
    type = 'Splicing'
  } else {
    type = 'TE'
  }
  all_logits_top$DataType = type
  return(all_logits_top)
}

all_logits_data = do.call(rbind, lapply(c(te_dist_files, ase_dist_files, splice_dist_files), read_dist_data))

si_data = filter(all_logits_data, Type != 'HallLabSV', Maf %in% c('novel')) %>%
  group_by(Riskratio,Lower,Upper,Pval,Type,Maf,Window,DataType) %>%
  sample_n(1) %>% ungroup()
sv_data = filter(all_logits_data, Type == 'HallLabSV', Maf %in% c('0-single')) %>%
  group_by(Riskratio,Lower,Upper,Pval,Type,Maf,Window,DataType) %>%
  sample_n(1) %>% ungroup()

si_data$Window = as.numeric(si_data$Window)
rp_cols = brewer.pal(9,'Blues')[c(4,6,8)]
names(rp_cols) = c('indels','SNPs', 'HallLabSV')
all.data = rbind(sv_data,si_data)
all.data$Method = factor(all.data$Method, levels=unique(all.data$Method))
all.data$Type = sapply(as.character(all.data$Type), function(x) 
  ifelse(x == 'HallLabSV', 'SVs', x))
all.data$Type = sapply(as.character(all.data$Type), function(x) 
  ifelse(x == 'SNPs', 'SNVs', x))
all.data$Type = factor(all.data$Type, levels=c('SVs', 'indels', 'SNVs'))
all.data$Significant = sapply(all.data$Pval, function(x)
  ifelse(x < 0.05, '1', '0'))

all.data$Method = sapply(all.data$Method, function(x) ifelse(x == 'ASE', 'aseOutliers',
                                                             ifelse(x == 'Splicing', 'sOutliers', 'eOutliers')))
all.data$Method = factor(all.data$Method, levels=c('eOutliers', 'aseOutliers', 'sOutliers'))
all.data$Type = factor(all.data$Type, levels=c('SNVs', 'indels', 'SVs'))
sf = ggplot(all.data, aes(x=Window, y=Riskratio,Group=Method)) + 
  geom_point(size=3,aes(color=Method), position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width=0,position=position_dodge(width=0.5)) +
  theme_bw() + ylab('Relative risk') + xlab('Distance downstream from gene (kb)') + 
  geom_line(aes(group=Method),size=0.5) +
  geom_hline(yintercept=1,color='grey') +
  scale_color_manual(values=dcols) + 
  gtex_v8_figure_theme() +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        strip.text.x=element_text(size=8),
        legend.title=element_blank()) +
  facet_wrap(Type~.,scales='free', ncol=3) + theme(strip.background = element_blank())


ggsave(sf, file=paste0(data_dir, 'paper_figures/revisions/supp/sfig_downstream.png'), width=7.2, height=4, units="in")

