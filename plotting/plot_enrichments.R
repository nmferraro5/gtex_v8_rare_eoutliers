library(ggplot2)
library(data.table)
library(dplyr)
library(RColorBrewer)

data_dir = '/Users/nicoleferraro/durga_local/data/goats_data/v8_data/'
risks = fread(paste0(data_dir, 'gtexV8.all.data.types.snps.indels.relative.risk.10kb.gnomad.txt'))

## if SNPs+indels
risks$MafBin = sapply(risks$Maf, function(x)
  ifelse(x == 'MAFnovel', 'novel',
         ifelse(x == 'MAFac1', 'single',
                ifelse(x == 'MAFac2', 'double',
                       ifelse(x == 'MAFgnomad0-1only', 'rare', 'low frequency')))))
## if SVs
risks$MafBin = sapply(risks$Maf, function(x)
  ifelse(x == 'MAF0-single', 'single',
         ifelse(x == 'MAFsingle-double', 'double',
                ifelse(x == 'MAF0-1', 'rare',
                       ifelse(x == 'MAF1-5', 'low frequency')))))

risks$MafBin = factor(risks$MafBin, levels=c('novel', 'single', 'double', 'rare', 'low frequency'))
risks$MafBin = factor(risks$MafBin, levels=c('single', 'double', 'rare', 'low frequency'))

si_plot = ggplot(risks, aes(x=MafBin,y=Riskratio,Group=Method)) +
  geom_point(size=6, aes(shape=Method),position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width=0,position=position_dodge(width=0.5)) +
  theme_bw() + ylab('Relative risk') + xlab('') +
  ggtitle('') + guides(shape=F) + 
  geom_hline(yintercept=1, linetype="dashed", colour="darkgrey") +
  theme(axis.title = element_text(size=24),
        axis.text.y = element_text(size=22),
        axis.text.x=element_text(size=22,angle=45,hjust=1),
        legend.title=element_blank(),
        legend.text=element_text(size=20),
        title=element_text(size=24)) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

### outlier counts
load(paste0(data_dir, 'outlier.counts.all.methods.RData'))
pvals = c('nonOutlier', pvals)
all_nums = rbind(cbind(te_nums, 'TE',pvals),cbind(splice_nums, 'Splicing',pvals),cbind(ase_nums, 'ASE',pvals))
all_nums = as.data.frame(all_nums)
all_nums$te_nums = as.numeric(as.character(all_nums$te_nums))
all_nums$pvals = factor(all_nums$pvals, levels=pvals)
mcols = c('#B8C5D6', '#A39BA8', '#272D2D')
ggplot(filter(all_nums, pvals != 'nonOutlier'), aes(x=pvals,y=te_nums,Group=V2)) + 
  geom_bar(stat='identity',aes(fill=V2),position='dodge') + theme_bw() +
  scale_fill_manual(values=mcols) + guides(fill=F) +
  xlab('Outlier threshold') + ylab('# outliers') +
  theme(axis.title = element_text(size=24),
        axis.text = element_text(size=22),
        legend.title=element_blank(),
        legend.text=element_text(size=20)) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

## plot single tissue
pvals = c(0.05, 0.01, 0.001, 0.0001, 0.00001, 0.0000001)
pnames = c('0.05', '0.01', '1e-03', '1e-04', '1e-05', '1e-07')
exp_single = fread(paste0(data_dir, 'enrichments/all.single.tissue.expression.relative.risks.txt')) %>%
  mutate(DataType = 'TotalExpression')
exp_single$Bin = rep(pnames, each=3)
splice_single = fread(paste0(data_dir, 'enrichments/all.single.tissue.splicing.relative.risks.txt')) %>%
  mutate(DataType = 'Splicing')
splice_single$Bin = rep(pnames, each=3)
ase_single = fread(paste0(data_dir, 'enrichments/all.single.tissue.ase.relative.risks.txt')) %>%
  mutate(DataType = 'ASE')
ase_single$Bin = rep(pnames, each=3)

both_single = rbind(exp_single, splice_single, ase_single)
both_single$Type = factor(both_single$Type, levels=c('SNPs','indels','SVs'))
both_single$DataType = factor(both_single$DataType, levels=c('ASE', 'Splicing','TotalExpression'))
both_single$Bin = factor(both_single$Bin, levels=pnames)
tcols = c('#335145', '#1E352F', '#598381')
tcols = c('#335145', '#C2C1A5', '#8C5E58')
single_plot = ggplot(both_single, aes(x=Bin,y=Riskratio,color=Type,shape=DataType)) +
  geom_point(size=6, position=position_dodge(width=0.7)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width=0, position=position_dodge(width=0.7)) +
  theme_bw() + ylab('Relative risk') + xlab('Outlier threshold') +
  ggtitle('') + scale_color_manual(values=tcols) + 
  geom_hline(yintercept=1, linetype="dashed", colour="darkgrey") +
  theme(axis.title = element_text(size=24),
        axis.text.y = element_text(size=22),
        axis.text.x=element_text(size=22,angle=45,hjust=1),
        legend.title=element_blank(),
        legend.text=element_text(size=20),
        title=element_text(size=24)) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

### plot tissue-specific enhancer enrichments
meth.colors = brewer.pal(9, 'YlGnBu')[3:9]
names(meth.colors) = c('EM','KNN','MEAN','PMD','SOFT','STFZ','MEDZ')
print(meth.colors)
enh_risks = fread(paste0(data_dir, 'enhancers/KNN.FDR.enhancer.snpindels.500kb.risks.jun1.txt'))
enh_risks$ZT = factor(enh_risks$ZT, levels=c(2,3,4,5,6,7,8,9))
enh_risks = filter(enh_risks, ZT %in% c(3,5,7))
enh_risks$Type = sapply(enh_risks$Type, function(x)
  ifelse(x == 'Any', 'Unmatched', x))
colnames(enh_risks)[6] = 'Zscore'
ggplot(enh_risks, aes(x=Type,y=Riskratio,Group=Zscore)) +
  geom_point(size=6, position=position_dodge(width=0.5),aes(alpha=Zscore),color=meth.colors[2]) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width=0, position=position_dodge(width=0.5)) +
  theme_bw() + ylab('Relative risk') + xlab('Enhancer type') +
  ggtitle('') + scale_alpha_manual(values=c(0.5,0.75,1)) +
  geom_hline(yintercept=1, linetype="dashed", colour="darkgrey") +
  theme(axis.title = element_text(size=24),
        axis.text = element_text(size=22),
        legend.title=element_text(size=20),
        legend.text=element_text(size=20),
        title=element_text(size=24)) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))



