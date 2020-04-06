library(data.table)
library(dplyr)
library(ggplot2)

data_dir = '/users/nferraro/data/goats_data/v8_data/'

gtex_gwas = fread(paste0(data_dir, 'coloc/gwas/gtex.traits.beta.gwas.overlap.txt'))
gtex_gwas$Chr = sapply(gtex_gwas$variant, function(x) strsplit(x, ':')[[1]][1])
gtex_gwas$Pos = sapply(gtex_gwas$variant, function(x) strsplit(x, ':')[[1]][2])
gtex_gwas$VID = paste(paste0('chr',gtex_gwas$Chr), gtex_gwas$Pos, sep=':')

exp_ws = fread(paste0(data_dir, 'coloc/gwas/expression.watershed.ukbb.variants.txt'))
sp_ws = fread(paste0(data_dir, 'coloc/gwas/splicing.watershed.ukbb.variants.txt'))
ase_ws = fread(paste0(data_dir, 'coloc/gwas/ase.watershed.ukbb.variants.txt'))

exp_ws = inner_join(exp_ws, gtex_gwas %>% select(minor_AF,beta,pval,Trait,VID), by='VID')
sp_ws = inner_join(sp_ws, gtex_gwas %>% select(minor_AF,beta,pval,Trait,VID), by='VID')
ase_ws = inner_join(ase_ws, gtex_gwas %>% select(minor_AF,beta,pval,Trait,VID), by='VID')

exp_ws_sum = exp_ws %>% group_by(VID) %>%
  mutate(MaxBeta = max(abs(beta))) %>%
  mutate(MaxWS = max(total_expression_watershed_posterior)) %>%
  sample_n(1) %>% ungroup() %>%
  select(VID,MaxBeta,MaxWS) %>% mutate(Type = 'TE')

sp_ws_sum = sp_ws %>% group_by(VID) %>%
  mutate(MaxBeta = max(abs(beta))) %>%
  mutate(MaxWS = max(splicing_watershed_posterior)) %>%
  sample_n(1) %>% ungroup() %>%
  select(VID,MaxBeta,MaxWS) %>% mutate(Type = 'Splicing')

ase_ws_sum = ase_ws %>% group_by(VID) %>%
  mutate(MaxBeta = max(abs(beta))) %>%
  mutate(MaxWS = max(ase_watershed_posterior)) %>%
  sample_n(1) %>% ungroup() %>%
  select(VID,MaxBeta,MaxWS) %>% mutate(Type = 'ASE')

ws_sum = rbind(exp_ws_sum, sp_ws_sum, ase_ws_sum) %>%
  group_by(VID) %>% mutate(MaxBetaAll = max(abs(MaxBeta))) %>%
  mutate(MaxWSAll = max(MaxWS)) %>% sample_n(1) %>% ungroup()

ggplot(ws_sum, aes(x=MaxBetaAll,y=MaxWSAll)) +
  geom_point() + theme_bw() +
  xlab('Max abs(beta)') + ylab('Max Watershed Score') +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        strip.text.x=element_text(size=16))

ws_all = rbind(exp_ws %>% mutate(posterior = total_expression_watershed_posterior) %>% select(-total_expression_watershed_posterior) %>% mutate(Type = 'TE'),
               sp_ws %>% mutate(posterior = splicing_watershed_posterior) %>% select(-splicing_watershed_posterior) %>% mutate(Type = 'Splicing'),
               ase_ws %>% mutate(posterior = ase_watershed_posterior) %>% select(-ase_watershed_posterior) %>% mutate(Type = 'ASE'))


ws_all$Gene = sapply(ws_all$sample_names, function(x) strsplit(x, ':')[[1]][2])
ws_all$Ind = sapply(ws_all$sample_names, function(x) strsplit(x, ':')[[1]][1])
ws_all = ws_all %>% select(-sample_names)
ws_all = ws_all %>% 
  mutate(WS_Bin = ifelse(posterior > 0.75, 0.75,
                         ifelse(posterior > 0.5, 0.5,
                                ifelse(posterior > 0.25, 0.25, 0))))
ws_all = ws_all %>% group_by(VID,Gene,Ind,Type) %>% top_n(1,abs(beta)) %>% ungroup()

ws_all$WS_Bin = factor(ws_all$WS_Bin, levels=c(0,0.25,0.5,0.75))
ggplot(ws_all, aes(x=Type, y=abs(beta),Group=WS_Bin)) + 
  geom_boxplot(aes(fill=WS_Bin)) + theme_bw() +
  xlab('') + ylab('Max abs(beta)') +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.title=element_text(size=20),
        legend.text=element_text(size=18)) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
## Intersect with outliers
medz_outliers = fread(paste0(data_dir, 'gtexV8.outliers.v8ciseQTLs.removed.medz.Z2.txt'))

sp_outliers = fread(paste0(data_dir, 'splicing/gtexV8.outliers.v8ciseQTLs.removed.splicing.Z2.txt'))

ase_data = fread(paste0(data_dir, 'ASE/combined.ad.scores.in.MEDIAN_4_10_update.tsv'))
ase_outliers = melt(ase_data) %>% filter(value < 2*pnorm(-abs(2))) 
colnames(ase_outliers)[1] = 'Gene'
gene_names = fread(paste0(data_dir, 'gencode.v26.GRCh38.genes.bed')) %>%
  select(-V1,-V2,-V3) %>% mutate(Gene = V4)
gene_names$Gene = sapply(gene_names$Gene, function(x) strsplit(x, '[.]')[[1]][1])
ase_outliers = inner_join(ase_outliers, gene_names, by='Gene') %>%
  mutate(GID = paste(variable,V4, sep='_'))

all_outliers = rbind(medz_outliers %>% mutate(value = MedZ) %>% select(-MedZ) %>% mutate(Type = 'TE'),
                     sp_outliers %>%  mutate(Type = 'Splicing'),
                     ase_outliers %>% mutate(Ind = variable,Gene = V4) %>% select(Gene,Ind,value) %>% mutate(Type = 'ASE'))


ws_all = merge(ws_all, all_outliers, by=c('Gene','Ind','Type'), all.x=T)
medz_genes = filter(ws_all, Type == 'TE', !is.na(value))$Gene
sp_genes = filter(ws_all, Type == 'Splicing', !is.na(value))$Gene
ase_genes = filter(ws_all, Type == 'ASE', !is.na(value))$Gene
ws_outliers = rbind(filter(ws_all, Type == 'TE', Gene %in% medz_genes),
               filter(ws_all, Type == 'Splicing', Gene %in% sp_genes),
               filter(ws_all, Type == 'ASE', Gene %in% ase_genes))
ws_outliers$IsOutlier = sapply(ws_outliers$value, function(x) ifelse(is.na(x), 0, 1))
ws_outliers$IsOutlier = factor(ws_outliers$IsOutlier, levels=c(0,1))
ws_outliers = ws_outliers %>% group_by(Gene,Ind,Type,VID,Trait) %>% sample_n(1) %>% ungroup()
ggplot(ws_outliers, aes(x=Type,y=posterior,Group=IsOutlier)) +
  geom_boxplot(aes(fill=IsOutlier)) + theme_bw() +
  xlab('') + ylab('Watershed posterior') +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.title=element_text(size=20),
        legend.text=element_text(size=18)) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ws_max = ws_outliers %>% filter(minor_AF < 0.05) %>%
  group_by(Gene,Type,Ind) %>% top_n(1,abs(beta)) %>% ungroup()
ggplot(ws_max, aes(x=Type,y=abs(beta),Group=IsOutlier)) +
  geom_boxplot(aes(fill=IsOutlier)) + theme_bw() +
  xlab('') + ylab('Max abs(beta)') +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.title=element_text(size=20),
        legend.text=element_text(size=18)) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ws_outliers = ws_outliers %>%
  mutate(OutlierWS = ifelse(IsOutlier == 1 & posterior > 0.25, 1, 0)) %>%
  filter(minor_AF < 0.01)
ws_outliers$OutlierWS = factor(ws_outliers$OutlierWS, levels=c(0,1))
ggplot(ws_outliers, aes(x=Type,y=abs(beta),Group=IsOutlier)) +
  geom_boxplot(aes(fill=IsOutlier)) + theme_bw() +
  xlab('') + ylab('abs(beta)') +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.title=element_text(size=20),
        legend.text=element_text(size=18)) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

### look at standing height for subset of genes with any outlier
height_data = filter(ws_all, Trait == '50_irnt')

height_top = height_data %>% group_by(Type,Gene) %>% top_n(1,abs(beta)) %>% ungroup()
height_other = height_data %>% filter(!(VID %in% height_top$VID))
height_plot = rbind(height_top %>% mutate(Cat='Top'), height_other %>% mutate(Cat = 'Other'))
ggplot(height_plot, aes(x=Type,y=posterior+0.00001)) + geom_boxplot(aes(fill=Cat)) + 
  theme_bw() + scale_y_log10() +
  ylab('Watershed posterior') + xlab('') +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.title=element_blank(),
        legend.text=element_text(size=18)) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

wilcox.test(filter(height_plot, Type == 'Splicing', Cat == 'Top')$posterior,
            filter(height_plot, Type == 'Splicing', Cat != 'Top')$posterior,
            alternative='g')


# use 2587 for top 10%
ws_high7 = ws_all %>% filter(posterior > 0.7) %>% filter(RankedBeta < 1294) %>%
  group_by(Trait) %>% mutate(NV = length(unique(VID))) %>%
  sample_n(1) %>% ungroup() %>% mutate(PT = 0.7)
ws_high8 = ws_all %>% filter(posterior > 0.8) %>% filter(RankedBeta < 1294) %>%
  group_by(Trait) %>% mutate(NV = length(unique(VID))) %>%
  sample_n(1) %>% ungroup() %>% mutate(PT = 0.8)
ws_high9 = ws_all %>% filter(posterior > 0.9) %>% filter(RankedBeta < 1294) %>%
  group_by(Trait) %>% mutate(NV = length(unique(VID))) %>%
  sample_n(1) %>% ungroup() %>% mutate(PT = 0.9)
ws_high = rbind(ws_high7, ws_high8, ws_high9)
ws_high$PT = factor(ws_high$PT, levels=c(0.7,0.8,0.9))
ws_high10 = ws_high %>% mutate(EBin = 'Top 10%')
ws_high5 = ws_high %>% mutate(EBin = 'Top 5%')
ws_high = rbind(ws_high10, ws_high5)
ggplot(ws_high, aes(x=PT,y=NV)) + geom_boxplot(aes(fill=EBin)) + theme_bw() +
  xlab('Watershed posterior threshold') + ylab('Number variants per trait') +
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=20),
        legend.text=element_text(size=20),
        legend.title=element_blank(),
        legend.key.size = unit(1.5, 'lines')) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ws_sum = ws_all %>% group_by(VID,Trait) %>% top_n(1,posterior) %>% ungroup() %>%
  group_by(VID) %>% top_n(1,-RankedBeta) %>% sample_n(1) %>% ungroup()

ggplot(ws_sum, aes(x=posterior, y=1/RankedBeta)) + geom_point() + theme_bw() +
  xlab('Watershed posterior') +
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=20)) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

## Take top variants per trait
ws_all$Gene = sapply(ws_all$sample_names, function(x) strsplit(x, ':')[[1]][2])
get_trait_data <- function(gtrait) {
  trait_data = filter(ws_all, Trait == gtrait)
  trait_top = trait_data %>% group_by(Type,Gene) %>% top_n(1,abs(beta)) %>% ungroup()
  trait_other = trait_data %>% filter(!(VID %in% trait_top$VID))
  trait_plot = rbind(trait_top %>% mutate(Cat='Top'), trait_other %>% mutate(Cat = 'Other'))
}

all_trait_compare = do.call(rbind, lapply(unique(ws_all$Trait), get_trait_data))
ggplot(all_trait_compare, aes(x=Type,y=posterior+0.00001)) + geom_boxplot(aes(fill=Cat)) + 
  theme_bw() + scale_y_log10() +
  ylab('Watershed posterior') + xlab('') +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.title=element_blank(),
        legend.text=element_text(size=18)) +
  annotate("text", x=1, y=1e-02, label='***',cex=5) +
  annotate("text", x=2, y=1e-02, label='***',cex=5) +
  annotate("text", x=3, y=1e-02, label='***',cex=5) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

wilcox.test(filter(all_trait_compare, Type == 'TE', Cat == 'Top')$posterior,
            filter(all_trait_compare, Type == 'TE', Cat != 'Top')$posterior,
            alternative='g')

### interesting examples
ws_all$Gene = sapply(ws_all$sample_names, function(x) strsplit(x, ':')[[1]][2])
ws_all = ws_all %>% group_by(Trait,Type) %>% mutate(RankedBeta = rank(-abs(beta))) %>% ungroup()
ws_top = filter(ws_all, posterior > 0.9, RankedBeta < 2587)
trait.table = fread(paste0(data_dir, 'coloc/gwas/ukbb.trait.table.txt'))
colnames(trait.table)[1] = 'Trait'
ws_top = merge(ws_top, trait.table, by='Trait')

cgene = 'ENSG00000107077'
gstart = 6757641
gend = 7175648
all_height_data = fread(paste0(data_dir, 'coloc/gwas/diabetes1.chr9.ukbb.sum.stats.tsv'))
all_height_data$Chr = sapply(all_height_data$variant, function(x) strsplit(x, ':')[[1]][1])
all_height_data$Pos = sapply(all_height_data$variant, function(x) strsplit(x, ':')[[1]][2])
all_height_data$VID = paste(paste0('chr', all_height_data$Chr), all_height_data$Pos, sep=':')
all_height_data = filter(all_height_data, Chr == 9, Pos > gstart - 10000, Pos < gend + 10000)
# eqtl for PUM3: chr9:2825014
# eqtl for TMEM18: chr2:678591
all_height_data$VOI = sapply(all_height_data$VID, function(x) ifelse(x == 'chr9:6916222', 'Outlier variant', 
                                                                     ifelse(x == 'chr9:7174433', 'eQTL', 'Other')))
all_height_data$VOI = factor(all_height_data$VOI, levels=c('Outlier variant','eQTL','Other'))
#all_height_data = all_height_data %>% filter(abs(beta) < 0.1)
ggplot(all_height_data, aes(x=minor_AF,y=abs(beta))) + 
  geom_point(aes(color=VOI,size=VOI)) + theme_bw() + guides(size=F) +
  geom_point(data = subset(all_height_data, VOI != 'Other'),
             aes(color = VOI, size = VOI)) +
  ggtitle('PIDD1 in Schizophrenia, chr11:802270 with Watershed posterior = 0.92') +
  scale_color_manual(values=c('hotpink', 'aquamarine', 'black')) + 
  scale_size_manual(values=c(4,4,2)) + 
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.title=element_blank(),
        legend.text=element_text(size=18),
        title=element_text(size=20)) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ex_data = filter(ws_all, GeneID == cgene, Description == 'Standing height') %>%
  group_by(VID) %>% top_n(1,posterior) %>% sample_n(1) %>% ungroup() %>%
  mutate(HighWS = ifelse(posterior > 0.9, 1, 0))
ex_data$HighWS = factor(ex_data$HighWS, levels=c(0,1))
ggplot(ex_data %>% filter(low_confidence_variant == 'FALSE'), aes(x=minor_AF,y=abs(beta))) + geom_point(aes(color=HighWS)) 




