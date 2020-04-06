### analysis for supplemental figure 2 in section 1 for TE - comparison of outlier types

library(gplots)
library(ggplot2)
library(ggthemes)
library(ggridges)
library(data.table)
library(dplyr)
library(RColorBrewer)
library(cowplot)
library(enrichR)

gtex_v8_figure_theme <- function() {
  return(theme(plot.title = element_text(face="plain",size=8), text = element_text(size=8),
               axis.text=element_text(size=7), panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(),panel.background = element_blank(), 
               axis.line = element_line(colour = "black"), 
               panel.border=element_blank(),
               legend.text = element_text(size=7), legend.title = element_text(size=8)))
}

data_dir = '/Users/nicoleferraro/durga_local/data/goats_data/v8_data/'

medz_data = fread(paste0(data_dir, 'gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.medz.filtered.txt'))
#medz_data = fread(paste0(data_dir, 'gtexV8.outliers.v8ciseQTLs.globalOutliers.removed.medz.txt'))
splice_data = melt(fread(paste0(data_dir, 'splicing/cross_tissue_covariate_method_none_no_global_outliers_ea_only_emperical_pvalue_gene_level.txt')))
ase_data = melt(fread(paste0(data_dir, 'ASE/median.ad.scores.uncorrected.no.global.outliers.v8.vg.tsv')))
colnames(ase_data) = c('GeneID', 'SampleName', 'DOT.score')
gene.mapper = fread(paste0(data_dir, 'gencode.v26.GRCh38.genes.bed'), header=F)
gene.mapper$GeneID = sapply(gene.mapper$V4, function(x) strsplit(x, '[.]')[[1]][1])
colnames(gene.mapper)[4] = 'Gene'
ase_data = inner_join(ase_data, gene.mapper %>% select(Gene,GeneID), by='GeneID')

colnames(splice_data)[1:2] = c('Gene', 'Ind')
colnames(ase_data)[2] = 'Ind'
combined_data = inner_join(medz_data, splice_data, by=c('Gene', 'Ind'))
combined_data = inner_join(combined_data, ase_data, by=c('Gene', 'Ind'))
combined_data = filter(combined_data, !is.na(DOT.score), !is.na(value))
shared_outliers = filter(combined_data, abs(MedZ) > 3, DOT.score < 0.0027, value < 0.0027)
symbols = fread(paste0(data_dir, 'gencode.v26.GRCh38.genes.symbols.txt'))
shared_symbols = filter(symbols, Gene %in% shared_outliers$Gene)
outlier_symbols = filter(symbols, Gene %in% combined_data$Gene)
exac = fread('/Users/nicoleferraro/durga_local/data/rds_data/exac_gene_metrics.txt')
exac_shared = filter(exac, gene %in% shared_symbols$Symbol)
exac_all = filter(exac, gene %in% outlier_symbols$Symbol)

te_sharing = c(nrow(filter(combined_data, abs(MedZ) > 3, DOT.score < 0.0027)) / nrow(filter(combined_data, abs(MedZ) > 3)), 
               nrow(filter(combined_data, abs(MedZ) > 3, value < 0.0027)) / nrow(filter(combined_data, abs(MedZ) > 3)))
ase_sharing = c(nrow(filter(combined_data, DOT.score < 0.0027, value < 0.0027)) / nrow(filter(combined_data, DOT.score < 0.0027)), 
               nrow(filter(combined_data, DOT.score < 0.0027, abs(MedZ) > 3)) / nrow(filter(combined_data, DOT.score < 0.0027)))
splice_sharing = c(nrow(filter(combined_data, value < 0.0027, DOT.score < 0.0027)) / nrow(filter(combined_data, value < 0.0027)), 
                nrow(filter(combined_data, value < 0.0027, abs(MedZ) > 3)) / nrow(filter(combined_data, value < 0.0027)))                  

all_sharing = as.data.frame(rbind(c(1, te_sharing), c(ase_sharing[1], 1, ase_sharing[2]), c(splice_sharing,1)))

sharing_df = melt(all_sharing)
sharing_df$Discovery = rep(c('eOutliers', 'aseOutliers', 'sOutliers'), 3)
sharing_df$Test = c('eOutliers', 'sOutliers', 'aseOutliers', 'aseOutliers', 'aseOutliers', 'eOutliers', 'sOutliers', 'eOutliers', 'sOutliers')

sharing_df$Test = factor(sharing_df$Test, levels=c('eOutliers', 'aseOutliers', 'sOutliers'))
sharing_df$Discovery = factor(sharing_df$Discovery, levels=c('eOutliers', 'aseOutliers', 'sOutliers'))

sfig_2A = ggplot(sharing_df, aes(y = Test, x = Discovery, fill = value)) +
  geom_tile(color = "grey") +
  xlab('Discovery method') + ylab('Test method') +
  geom_text(aes(label = round(value, 2)),cex=3) +
  scale_fill_gradient("Sharing \nfraction",
                      low = "white", high = "#BC6F8B",
                      na.value = "black", limits = c(0, .52)) +
  gtex_v8_figure_theme()

combined_data$Cat = sapply(1:nrow(combined_data), function(x) 
  ifelse(combined_data$DOT.score[x] < 0.0027 && abs(combined_data$MedZ[x]) > 3, 'Shared Outlier (n=319)', 
         ifelse(combined_data$DOT.score[x] < 0.0027, 'aseOutlier (n=1890)',
                ifelse(abs(combined_data$MedZ[x]) > 3, 'eOutlier (n=305)', 'Non-outlier (n=496536)'))))

gkeep = unique(filter(combined_data, Cat != 'Non-outlier (n=496536)')$Gene)
sharing_data = filter(combined_data, Gene %in% gkeep)
sharing_data$Cat = factor(sharing_data$Cat, levels=c('Shared Outlier (n=319)', 'eOutlier (n=305)', 'aseOutlier (n=1890)', 'Non-outlier (n=496536)'))
sfig_2C = ggplot(sharing_data, aes(y = Cat, x = abs(MedZ))) + 
  stat_density_ridges(quantiles=2, quantile_lines=T) +
  gtex_v8_figure_theme() + ylab('') + xlab('|Median Z-score|') 

### sharing annotations
ase_annot = fread(paste0(data_dir, 'ASE/gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.ase.variantAnnotations.v8.vg.txt'))
ase_annot$Ind = paste0('GTEX-', ase_annot$indiv_id)
ase_only = filter(sharing_data, Cat == 'aseOutlier (n=1890)')
colnames(ase_annot)[2] = 'Gene'
ase_only = inner_join(ase_only, ase_annot, by=c('Ind', 'Gene'))
ase_only$medz_bin = sapply(ase_only$MedZ, function(x) ifelse(x > 2, '2~3',
                                                             ifelse(x > 1, '1~2',
                                                                    ifelse(x > -1, '-1~1',
                                                                           ifelse(x > -2, '-2~-1', '-3~-2')))))
ase_summary = ase_only %>% group_by(medz_bin) %>% mutate(NTotal = n()) %>% ungroup() %>%
  group_by(medz_bin,variant_cat) %>% mutate(NBin = n()) %>% mutate(NProp = NBin / NTotal) %>%
  sample_n(1) %>% ungroup() %>% filter(!is.na(variant_cat)) %>%
  mutate(CatCol = rgb(variant_color1, variant_color2, variant_color3))

plot_cols = unique(ase_summary$CatCol)
names(plot_cols) = unique(ase_summary$variant_cat)

ase_summary$medz_bin = factor(ase_summary$medz_bin,
                               levels=c('-3~-2','-2~-1','-1~1','1~2','2~3'))
ase_summary$variant_cat = factor(ase_summary$variant_cat,
                                  levels=c('no_variant','other_noncoding','coding','conserved_noncoding','TSS','stop','frameshift','splice','TE','INV','BND','DEL','CNV','DUP'))

sfig_2D = ggplot(ase_summary, aes(x=medz_bin,y=NProp,Group=variant_cat)) +
  geom_bar(aes(fill=variant_cat),color='black',stat='identity') +
  theme_bw() + coord_flip() + ylab('') + xlab('') +
  scale_fill_manual(values=plot_cols) +
  ggtitle('aseOutliers') + xlab('Median Z-score bin') +
  gtex_v8_figure_theme() +
  theme(legend.title=element_blank(),
        legend.key.height = unit(0.1, "cm"),
        legend.key.width = unit(0.1, "cm"))

shared_annot = as.data.frame(table(inner_join(shared_outliers, ase_annot, by=c('Ind', 'Gene'))$variant_cat))
shared_annot$Var1 = factor(shared_annot$Var1, levels=c('no_variant', 'DEL', 'conserved_noncoding', 'other_noncoding', 'splice'))
sfig_2B = ggplot(shared_annot, aes(x=Var1,y=Freq/35)) + geom_bar(stat='identity') +
  xlab('') + ylab('Proportion with nearby rare variant') +
  gtex_v8_figure_theme() + 
  theme(axis.text.x=element_text(hjust=1,angle=45))

first_row <- plot_grid(sfig_2A, sfig_2B, labels = c('A','B'), ncol=2, align='h')
second_row <- plot_grid(sfig_2C, sfig_2D, labels = c('C','D'), ncol=2, align='h')
sfig_2 = plot_grid(first_row, second_row, nrow = 2, align='v')

ggsave(sfig_2, file=paste0(data_dir, 'paper_figures/sfig_sharing_ase.pdf'), width=7.2, height=7.2,units="in")

## extremes table
mpercentiles = quantile(abs(filter(medz_data, abs(MedZ) > 3)$MedZ), probs=seq(0.9,1,0.01))
apercentiles = quantile(-log10(filter(ase_data, DOT.score < 0.0027)$DOT.score), probs=seq(0.9,1,0.01))
smin = min(filter(splice_data, value != 0)$value)
splice_data$adjvalue = sapply(splice_data$value, function(x) ifelse(x == 0, 5e-08, x))
spercentiles = quantile(-log10(filter(splice_data, adjvalue < 0.0027)$adjvalue), probs=seq(0.9,1,0.01))

medz_extremes = filter(medz_data, abs(MedZ) > mpercentiles[10])
ase_extremes = filter(ase_data, -log10(DOT.score) > apercentiles[10])
splice_extremes = filter(splice_data, -log10(adjvalue) >= spercentiles[10])
extreme_genes = unique(c(medz_extremes$Gene, ase_extremes$Gene, splice_extremes$Gene))
gencode = read.table(paste0(data_dir, 'gencode.v26.GRCh38.genes.gtf'),fill=NA) %>%
  select(V13,V19)
extreme_symbols = unique(filter(gencode, V13 %in% extreme_genes)$V19)






