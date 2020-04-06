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

## high cholesterol
cgene = 'ENSG00000075234.16'
gstart = 46267961
gend = 46294008

all_trait_data = fread(paste0(data_dir, 'coloc/gwas/traits/cholesterol.chr22.ukbb.sum.stats.tsv'))
all_trait_data$Chr = sapply(all_trait_data$variant, function(x) strsplit(x, ':')[[1]][1])
all_trait_data$Pos = sapply(all_trait_data$variant, function(x) strsplit(x, ':')[[1]][2])
all_trait_data$VID = paste(paste0('chr', all_trait_data$Chr), all_trait_data$Pos, sep=':')

all_trait_data$Pos = as.numeric(all_trait_data$Pos)
all_trait_plot_data = filter(all_trait_data, Chr == 22, Pos > gstart - 10000, Pos < gend + 10000)
all_trait_plot_data$VOI = sapply(all_trait_plot_data$VID, function(x) ifelse(x == 'chr22:46291947', 'Watershed variant', 'Other'))
all_trait_plot_data$VOI = factor(all_trait_plot_data$VOI, levels=c('Watershed variant','Other'))
cc_ratio = 43957/317184
scale_factor = cc_ratio * (1-cc_ratio)
all_trait_plot_data$beta_scaled = all_trait_plot_data$beta / scale_factor

sfigA = ggplot(all_trait_plot_data %>% filter(low_confidence_variant == 'FALSE'), aes(x=minor_AF,y=abs(beta_scaled))) +
  geom_point(aes(color=VOI),size=1) + theme_bw() + guides(size=F) +
  geom_point(data = subset(all_trait_plot_data, VOI != 'Other'),
             aes(color = VOI), size=2) +
  ylab('|Effect size|') +
  scale_color_manual(values=c('#CC00B4', 'black')) + xlab('UKBB MAF') +
  gtex_v8_figure_theme() + guides(color=F)

### MVP replication
tg_data = fread(paste0(data_dir, 'coloc/gwas/MVP/CLEANED.MVP.Euro.TG.results.EQC.csv.gz_rs564796245_250kb.txt')) %>%
  mutate(Trait = 'Triglycerides')
hdl_data = fread(paste0(data_dir, 'coloc/gwas/MVP/CLEANED.MVP.Euro.HDL.results.EQC.csv.gz_rs564796245_250kb.txt')) %>%
  mutate(Trait = 'HDL')
ldl_data = fread(paste0(data_dir, 'coloc/gwas/MVP/CLEANED.MVP.Euro.LDL.results.EQC.csv.gz_rs564796245_250kb.txt')) %>%
  mutate(Trait = 'LDL')
tc_data = fread(paste0(data_dir, 'coloc/gwas/MVP/CLEANED.MVP.Euro.TC.results.EQC.csv.gz_rs564796245_250kb.txt')) %>%
  mutate(Trait = 'Total Cholesterol')

all_chol_data = rbind(tg_data, hdl_data, ldl_data, tc_data)
all_chol_data$VID = paste(all_chol_data$CHR, all_chol_data$BP, sep=':')
voi = '22:46687844'
all_chol_data$VOI = sapply(all_chol_data$VID, function(x) ifelse(x == voi, 1, 0))

all_chol_data = fread(paste0(data_dir, 'coloc/gwas/MVP/MVP.chol.results.txt'))
all_chol_data$AF_NFE = as.numeric(all_chol_data$AF_NFE)
all_chol_data$GAF = sapply(all_chol_data$AF_NFE, function(x) ifelse(is.na(x), 0, 
                                                                    ifelse(x > 0.5, 1-x, x)))
all_chol_data$VOI = factor(all_chol_data$VOI, levels=c(0,1))
all_chol_data$GAF = as.numeric(all_chol_data$GAF)
dat_text <- data.frame(
  Trait = c("HDL", "LDL", "Total Cholesterol", "Triglycerides"),
  label   = c('99%', '94%', '97%', '94%'))
chol_rare_data = all_chol_data %>% filter(GAF < 0.001) %>% group_by(Trait) %>% 
  mutate(VPercentile = ntile(abs(BETA), 100)) %>% ungroup()
sfigB = ggplot(all_chol_data, aes(x=GAF,y=abs(BETA))) +
  geom_point(aes(color=VOI),size=0.4) + theme_bw() + 
  geom_point(data = subset(all_chol_data, VOI == 1),
             aes(color = VOI),size=2) +
  scale_color_manual(values=c('black', '#CC00B4')) +
  xlab('gnomAD minor allele frequency') + guides(color=F,size=F) +
  gtex_v8_figure_theme() + facet_wrap('Trait', ncol=2) +
  ylab('|Effect size|') + 
  theme(strip.background = element_blank(),
        strip.text.x=element_text(size=8))

fig_dir = '/Users/nicoleferraro/durga_local/data/goats_data/v8_data/paper_figures/'
chol_fig <- plot_grid(sfigA, sfigB, labels = c('A', 'B'), ncol=2, axis='tblr', align='h', rel_widths=c(1,2))
ggsave(sfigB, file=paste0(fig_dir, 'sfig_topmed_revisions.pdf'), width=7.2, height=4, units="in")

### MVP replication - BMI
bmi_data_v1 = fread(paste0(data_dir, 'coloc/gwas/MVP/BMI-height_rare_variants_in_MVP_withCHR.POS/rs138799125_bmi.EUR.txt')) %>%
  mutate(Trait = 'BMI', Var = 'rs138799125')
bmi_data_v2 = fread(paste0(data_dir, 'coloc/gwas/MVP/BMI-height_rare_variants_in_MVP_withCHR.POS/rs149737009_bmi.EUR.txt')) %>%
  mutate(Trait = 'BMI', Var = 'rs149737009')
bmi_data_v3 = fread(paste0(data_dir, 'coloc/gwas/MVP/BMI-height_rare_variants_in_MVP_withCHR.POS/rs184629044_bmi.EUR.txt')) %>%
  mutate(Trait = 'BMI', Var = 'rs184629044')
bmi_data_v4 = fread(paste0(data_dir, 'coloc/gwas/MVP/BMI-height_rare_variants_in_MVP_withCHR.POS/rs192797357_bmi.EUR.txt')) %>%
  mutate(Trait = 'BMI', Var = 'rs192797357')
bmi_data_v5 = fread(paste0(data_dir, 'coloc/gwas/MVP/BMI-height_rare_variants_in_MVP_withCHR.POS/rs35835984_bmi.EUR.txt')) %>%
  mutate(Trait = 'BMI', Var = 'rs35835984')

bmi_vars = fread(paste0(data_dir, 'paper_tables/high.watershed.effect.size.variants.old.txt')) %>% filter(grepl('BMI', Description))
bmi_vars = filter(bmi_vars, variant %in% c('11:840224:T:C', '11:825341:C:T', '11:846407:C:T'))
bmi.chr11.stats = fread(paste0(data_dir, 'coloc/gwas/traits/bmi.chr11.ukbb.sum.stats.tsv'))
bmi.chr11.stats$VOI = sapply(bmi.chr11.stats$variant, function(x) ifelse(x %in% bmi_vars$variant, 1, 0))
dlow = 825341 - 250000
dhigh = 846407 + 250000
bmi.chr11.stats$Pos = sapply(bmi.chr11.stats$variant, function(x) as.numeric(strsplit(x, ':')[[1]][2]))
bmi.chr11.plot = filter(bmi.chr11.stats, Pos > dlow, Pos < dhigh, low_confidence_variant == FALSE)
bmi.chr11.plot$VOI = factor(bmi.chr11.plot$VOI, levels=c(0,1))

ggplot(bmi.chr11.plot, aes(x=minor_AF,y=beta)) +
  geom_point(aes(color=VOI),size=0.75) + theme_bw() + 
  geom_point(data = subset(bmi.chr11.plot, VOI == 1),
             aes(color = VOI),size=2.5) +
  scale_color_manual(values=c('black', '#CC00B4')) +
  xlab('UKBB MAF') + guides(color=F,size=F) +
  ylab('|Effect size|') + 
  theme(strip.background = element_blank(),
        strip.text.x=element_text(size=8),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14))

bmi_data_v3 = rbind(bmi_data_v3, list(11, 825341, -0.0114851, 0.0268615, 'variant', 'BMI', 'rs184629044'))
bmi_data_v4 = rbind(bmi_data_v4, list(11, 846407, -0.0117553, 0.0269907, 'variant', 'BMI', 'rs192797357'))
bmi_data_v5 = rbind(bmi_data_v5, list(11, 840224, 0.00169545, 0.0295364, 'variant', 'BMI', 'rs35835984'))

gafs = fread(paste0(data_dir, 'coloc/gwas/MVP/11.hg19.all.INFO'))
colnames(gafs)[1] = 'CHR'

bmi_mvp_data = rbind(bmi_data_v3, bmi_data_v4, bmi_data_v5) %>% select(CHR,POS,BETA,SE,location) %>% distinct()
bmi_mvp_data = merge(bmi_mvp_data, gafs, by=c('CHR', 'POS'), all.x=T)

bmi_mvp_data$AF_NFE = sapply(bmi_mvp_data$AF_NFE, function(x) ifelse(is.na(x), 0, x))
bmi_mvp_data$AF_NFE = sapply(bmi_mvp_data$AF_NFE, function(x) strsplit(x, ',')[[1]][1])
bmi_mvp_data$AF_NFE = as.numeric(bmi_mvp_data$AF_NFE)
bmi_mvp_data$VOI = sapply(bmi_mvp_data$location, function(x) ifelse(x == 'variant', 1, 0))
bmi_mvp_data$VOI = factor(bmi_mvp_data$VOI, levels=c(0,1))
bmi_mvp_data$AF_NFE = sapply(bmi_mvp_data$AF_NFE, function(x) ifelse(x > 0.5, 1-x, x))

rare_bmi_mvp = filter(bmi_mvp_data, AF_NFE < 0.01) %>% mutate(VPercentile = ntile(abs(BETA), 100))

ggplot(bmi_mvp_data, aes(x=AF_NFE,y=abs(BETA))) +
  geom_point(aes(color=VOI),size=0.75) + theme_bw() + 
  geom_point(data = subset(bmi_mvp_data, VOI == 1),
             aes(color = VOI),size=2.5) +
  scale_color_manual(values=c('black', '#CC00B4')) +
  xlab('gnomAD MAF') + guides(color=F,size=F) +
  ylab('|Effect size|') + 
  theme(strip.background = element_blank(),
        strip.text.x=element_text(size=8),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14))

### look at outlier replication with MESA
rare_snps = fread(paste0(data_dir, 'coloc/gwas/MESA/gtex.outliers.mesa.zscores.snps.txt'))
rare_snps = rare_snps %>% select(Ind,Gene,TOPID,MedZ,value) %>% distinct()
rare_indels = fread(paste0(data_dir, 'coloc/gwas/MESA/gtex.outliers.mesa.zscores.indels.txt'))
rare_indels = rare_indels %>% select(Ind,Gene,TOPID,MedZ,value) %>% distinct()
prop_snp = nrow(filter(rare_snps, abs(value) > 2)) / nrow(rare_snps)
prop_indel = nrow(filter(rare_indels, abs(value) > 2)) / nrow(rare_indels)

ggplot(rare_snps, aes(x=MedZ, y=value)) + geom_point() + theme_bw() +
  geom_hline(yintercept=2, color='firebrick4') +
  geom_hline(yintercept=-2, color='firebrick4') +
  ggtitle('Rare SNPs') +
  xlab('MESA Z-score') + ylab('GTEx Z-score') +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        title=element_text(size=16))

ggplot(rare_indels, aes(x=MedZ, y=value)) + geom_point() + theme_bw() +
  geom_hline(yintercept=2, color='firebrick4') +
  geom_hline(yintercept=-2, color='firebrick4') +
  xlab('MESA Z-score') + ylab('GTEx Z-score') +
  ggtitle('Rare indels') +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        title=element_text(size=16))

## watershed plot
mean_values = fread(paste0(data_dir, 'coloc/gwas/MESA/TOPMed_MESA_RNAseq_Pilot_RNASeQCv2.0.0.gene.peer.ztrans.meanExams.txt'))

rare_ws_snps = fread(paste0(data_dir, 'coloc/gwas/MESA/gtex.outliers.mesa.zscores.snps.watershed.txt'))
rare_ws_snps = inner_join(rare_ws_snps, mean_values, by=c('Id', 'variable'))

ws_snps = rare_ws_snps %>% select(MedianWS,MeanValue,Id) %>% filter(MedianWS < 0.01 | MedianWS > 0.8) %>%
  mutate(Category = 'Watershed') %>%
  mutate(Bin = ifelse(MedianWS <= 0.01, 'Median GTEx posterior < 0.01', 'Median GTEx posterior > 0.8'))
gam_snps = rare_ws_snps %>% select(MedianGAM,MeanValue,Id) %>% filter(MedianGAM < 0.01 | MedianGAM > 0.8) %>%
  mutate(Category = 'GAM') %>%
  mutate(Bin = ifelse(MedianGAM <= 0.01, 'Median GTEx posterior < 0.01', 'Median GTEx posterior > 0.8'))
colnames(ws_snps)[1] = 'MedianPosterior'
colnames(gam_snps)[1] = 'MedianPosterior'
ws_keep = filter(ws_snps, Bin == 'Median GTEx posterior > 0.8') %>% select(Id)
gam_keep = filter(gam_snps, Bin == 'Median GTEx posterior > 0.8') %>% select(Id)

both_data = rbind(ws_snps %>% filter(Id %in% ws_keep$Id), gam_snps %>% filter(Id %in% gam_keep$Id)) %>%
  distinct() %>% group_by(MeanValue,Id,Category,Bin) %>%
  sample_n(1) %>% ungroup()

both_data$LogPval = -log10(2*pnorm(-abs(both_data$MeanValue)))


ggplot(both_data, aes(x=Bin,y=abs(MeanValue))) + geom_boxplot(aes(fill=Category)) + theme_bw() +
  scale_fill_manual(values=c('firebrick4', 'steelblue3')) +
  xlab('') + ylab('|Z-score|') + theme(legend.title=element_blank())











