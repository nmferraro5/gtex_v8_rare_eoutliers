## check variants in TOPMed

library(ggplot2)
library(data.table)
library(dplyr)
library(RColorBrewer)
library(cowplot)


gtex_v8_figure_theme <- function() {
  return(theme(plot.title = element_text(face="plain",size=8), text = element_text(size=8),
               axis.text=element_text(size=7), panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(),panel.background = element_blank(), 
               panel.border = element_blank(),
               axis.line = element_line(colour = "black"), 
               legend.text = element_text(size=7), legend.title = element_text(size=8)))
}

data_dir = '/Users/nicoleferraro/durga_local/data/goats_data/v8_data/'

cgene = 'ENSG00000075234.16'
gstart = 46267961
gend = 46294008

### MVP replication
tg_data = fread(paste0(data_dir, 'topmed/chr22_46687844.250kb.GTEXLookup.JHS.trigs_log.txt')) %>%
  mutate(Trait = 'Triglycerides')
hdl_data = fread(paste0(data_dir, 'topmed/chr22_46687844.250kb.GTEXLookup.JHS.hdl.txt')) %>%
  mutate(Trait = 'HDL')
ldl_data = fread(paste0(data_dir, 'topmed/chr22_46687844.250kb.GTEXLookup.JHS.ldl_adj.txt')) %>%
  mutate(Trait = 'LDL')
tc_data = fread(paste0(data_dir, 'topmed/chr22_46687844.250kb.GTEXLookup.JHS.totchol_adj.txt')) %>%
  mutate(Trait = 'Total Cholesterol')

vids = paste0('chr', paste(tg_data$CHROM, paste(tg_data$BEG, tg_data$END, sep='-'), sep=':'))
vids_hg38 = fread(paste0(data_dir, 'topmed/topmed_chol_variants_hg38.bed'), header=F)
vids_failed = fread(paste0(data_dir, 'topmed/variants_failed_liftover.txt'), header=F)

vid_map = cbind(vids_hg38, vids[which(!(vids %in%  vids_failed$V1))])
all_chol_data = rbind(tg_data, hdl_data, ldl_data, tc_data)
all_chol_data$HG19_VID = rep(vids, 4)
all_chol_data = filter(all_chol_data, HG19_VID %in% vid_map$V2)
colnames(vid_map) = c('HG38_VID', 'HG19_VID')
all_chol_data = inner_join(all_chol_data, vid_map, by='HG19_VID')
all_chol_data$HG38_Pos = sapply(all_chol_data$HG38_VID, function(x) as.numeric(strsplit(strsplit(x, ':')[[1]][2], '-')[[1]][1]))

gnomad_afs = fread(paste0(data_dir, 'topmed/22.all.INFO'))
gnomad_afs$CHROM = as.numeric(gnomad_afs$CHROM)
colnames(gnomad_afs)[2] = 'HG38_Pos'
all_chol_data = merge(all_chol_data, gnomad_afs, by=c('CHROM', 'HG38_Pos'), all.x=T)
all_chol_data$HG38_Pos = as.numeric(all_chol_data$HG38_Pos)
all_chol_data$AF_NFE = as.numeric(all_chol_data$AF_NFE)
all_chol_data$AF_NFE = sapply(all_chol_data$AF_NFE, function(x) ifelse(is.na(x), 0, x))

all_chol_data$VOI = sapply(all_chol_data$HG38_Pos, function(x) ifelse(x == 46291947, 1, 0))
all_chol_data$VOI = factor(all_chol_data$VOI, levels=c(0,1))
sfig = ggplot(all_chol_data %>% filter(AF_NFE <= 0.5), aes(x=AF_NFE,y=abs(BETA))) +
  geom_point(aes(color=VOI),size=0.4) + theme_bw() + 
  geom_point(data = subset(all_chol_data, VOI == 1),
             aes(color = VOI),size=2) +
  scale_color_manual(values=c('black', '#CC00B4')) +
  xlab('gnomAD minor allele frequency') + guides(color=F,size=F) +
  gtex_v8_figure_theme() + facet_wrap('Trait', ncol=2) +
  ylab('|Effect size|') + 
  theme(strip.background = element_blank(),
        strip.text.x=element_text(size=8))

## filtering to UKBB variants
all_trait_data = fread(paste0(data_dir, 'coloc/gwas/traits/cholesterol.chr22.ukbb.sum.stats.tsv'))
all_trait_data$Pos = sapply(all_trait_data$variant, function(x) as.numeric(strsplit(x, ':')[[1]][2]))
chol_ukbb_data = filter(all_chol_data, HG38_Pos %in% all_trait_data$Pos)
sfig2 = ggplot(chol_ukbb_data %>% filter(AF_NFE <= 0.5), aes(x=AF_NFE,y=abs(BETA))) +
  geom_point(aes(color=VOI),size=0.4) + theme_bw() + 
  geom_point(data = subset(all_chol_data, VOI == 1),
             aes(color = VOI),size=2) +
  scale_color_manual(values=c('black', '#CC00B4')) +
  xlab('gnomAD minor allele frequency') + guides(color=F,size=F) +
  gtex_v8_figure_theme() + facet_wrap('Trait', ncol=2) +
  ylab('|Effect size|') + 
  theme(strip.background = element_blank(),
        strip.text.x=element_text(size=8))

ggsave(sfig, file=paste0(data_dir, 'paper_figures/sfig_jhs_chol.pdf'), width=7.2, height=7.2, units="in")


chol_rare_data = all_chol_data %>% filter(AF_NFE < 0.001) %>% group_by(Trait) %>% 
  mutate(VPercentile = ntile(abs(BETA), 100)) %>% ungroup()

chol_rare_uk_data = chol_ukbb_data %>% filter(AF_NFE < 0.01) %>% group_by(Trait) %>% 
  mutate(VPercentile = ntile(abs(BETA), 100)) %>% ungroup()


