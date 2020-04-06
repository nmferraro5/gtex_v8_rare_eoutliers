library(data.table)
library(dplyr)
library(ggplot2)
library(cowplot)

data_dir = '/Users/nicoleferraro/Documents/Stanford/MontgomeryLab/mpra-v2-master/output/'
gtex_dir = '/Users/nicoleferraro/durga_local/data/goats_data/v8_data/'

gtex_v8_figure_theme <- function() {
  return(theme(plot.title = element_text(face="plain",size=8), text = element_text(size=8),
               axis.text=element_text(size=7), panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(),panel.background = element_blank(), 
               panel.border = element_blank(),
               axis.line = element_line(colour = "black"), 
               legend.text = element_text(size=7), legend.title = element_text(size=8)))
}

deseq2_results = fread(paste0(gtex_dir, 'mpra/GTEX_DESeq2_results.txt'))
deseq2_results$Gene = sapply(deseq2_results$geneID, function(x) strsplit(x, '=')[[1]][2])
deseq2_results$Chr = sapply(deseq2_results$chr, function(x) strsplit(x, '=')[[1]][2])
deseq2_results$Pos = sapply(deseq2_results$pos, function(x) strsplit(x, '=')[[1]][2])

deseq2_results = deseq2_results %>% select(Gene,Chr,Pos,baseMean_allele,log2FoldChange_allele,pvalue_allele,padj_allele)
deseq2_results$VID = paste(deseq2_results$Chr, deseq2_results$Pos, sep=':')
deseq2_results = deseq2_results %>% group_by(Gene,Chr,Pos) %>% top_n(1,-padj_allele) %>% sample_n(1) %>% ungroup()

ws_scores = fread(paste0(gtex_dir, 'mpra/mpra_deseq_watershed_genes.txt'))
ws_converted = fread(paste0(gtex_dir, 'mpra/watershed_variants_converted_positions.txt'))
colnames(ws_converted)[1:2] = c('Chr', 'Pos')
ws_scores = inner_join(ws_scores, ws_converted, by=c('Chr','Pos'))
ws_scores$VID = paste(ws_scores$CHR_19, ws_scores$Pos_19, sep=':')
ws_summary = ws_scores %>% group_by(Gene,VID) %>% top_n(1,total_expression_watershed_posterior) %>% 
  sample_n(1) %>% ungroup()
colnames(deseq2_results)[1] = 'GeneID'
deseq_watershed = inner_join(deseq2_results, ws_summary, by=c('VID', 'GeneID'))
deseq_watershed$HighWS = sapply(deseq_watershed$total_expression_watershed_posterior, function(x) ifelse(x > 0.5, 'High', 'Low'))
ggplot(deseq_watershed, aes(x=total_expression_watershed_posterior, y=abs(log2FoldChange_allele))) + 
  geom_point() + theme_bw() + theme(axis.text=element_text(size=12),
                                    axis.title=element_text(size=14))
ggplot(deseq_watershed, aes(x=HighWS, y=abs(log2FoldChange_allele))) + 
  geom_boxplot() + theme_bw() + theme(axis.text=element_text(size=12),
                                    axis.title=element_text(size=14)) +
  xlab('Watershed posterior bin') +
  annotate("text", x=1.5, y=1.5, label='p=0.07',cex=5)

outliers = fread(paste0(gtex_dir, 'gtexV8.outliers.v8ciseQTLs.globalOutliers.removed.medz.txt'))
deseq_watershed_outliers = filter(deseq_watershed, Gene %in% outliers$Gene)
top_genes = filter(deseq_watershed_outliers, total_expression_watershed_posterior > 0.5)
ggplot(deseq_watershed_outliers, aes(x=total_expression_watershed_posterior, y=abs(log2FoldChange_allele))) + 
  geom_point() + theme_bw() + theme(axis.text=element_text(size=12),
                                    axis.title=element_text(size=14))
deseq_watershed_outliers$WS_Bin = sapply(deseq_watershed_outliers$total_expression_watershed_posterior, 
                                         function(x) ifelse(x >= 0.5, 'High','Low'))
sf1 = ggplot(deseq_watershed_outliers, aes(x=WS_Bin, y=abs(log2FoldChange_allele))) + 
  geom_boxplot() + theme_bw() + 
  gtex_v8_figure_theme() + scale_y_log10() +
  xlab('Watershed posterior bin') + ylab('|log(allelic fold-change)|') +
  annotate("text", x=1.5, y=1, label='p=0.026',cex=3)

sf2 = ggplot(deseq_watershed_outliers, aes(x=WS_Bin, y=-log10(padj_allele))) + 
  geom_boxplot() + theme_bw() + 
  gtex_v8_figure_theme() + scale_y_log10() +
  xlab('Watershed score') + ylab('-log10(p adjusted)') +
  annotate("text", x=1.5, y=1, label='p=0.084',cex=3)

fig_1 = plot_grid(sf1, sf2, ncol = 2, labels=c('A', 'B'), align='hv', axis='tlbr')
ggsave(sf1, file=paste0(data_dir, 'paper_figures/revisions/supp/sfig_mpra.pdf'), width=7.2, height=5, units="in")

medz_outliers = fread(paste0(gtex_dir, 'gtexV8.outliers.v8ciseQTLs.removed.medz.Z2.withrarevars.txt'))
medz_outliers$Gene = sapply(medz_outliers$Gene, function(x) strsplit(x, '[.]')[[1]][1])
medz_outliers$VID = paste(medz_outliers$Chr, medz_outliers$Pos, sep=':')
converted_vars = fread(paste0(data_dir, 'hglft_genome_24b0b_4bb8e0.bed'))

medz_outliers = filter(medz_outliers, Gene %in% deseq2_results$Gene)
medz_outliers = cbind(medz_outliers, converted_vars)
medz_outliers$VID = paste(medz_outliers$V1,medz_outliers$V2,sep=':')

deseq2_results$IsOutlier = sapply(deseq2_results$VID, function(x) ifelse(x %in% medz_outliers$VID, 1, 0))
deseq2_outliers = filter(deseq2_results, GeneID %in% medz_outliers$Gene)
deseq2_outliers$IsOutlier = factor(deseq2_outliers$IsOutlier, levels=c(0,1))
deseq2_outliers$IsSignificant = sapply(deseq2_outliers$padj_allele, function(x) ifelse(x < 0.05, 1, 0))
deseq2_outliers$IsSignificant = factor(deseq2_outliers$IsSignificant, levels=c(0,1))

ggplot(deseq2_outliers, aes(x=baseMean_allele,y=log2FoldChange_allele)) + 
  geom_point(aes(shape=IsOutlier,color=IsSignificant),size=3) +
  scale_color_manual(values=c('black','red')) +
  xlab('Mean expression') + ylab('Log2(Allelic fold change)') + theme_bw() +
  annotate("text", x=55000,y=2.1,label='COA3, chr17:42798769, MedZ = -7.13, Watershed = 0.94',cex=5) +
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=16),
        legend.text=element_text(size=16),
        legend.title=element_text(size=16),
        legend.position=c(0.9,0.25)) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

colnames(deseq2_results)[1] = 'Gene'
deseq2_results$Pos = as.numeric(deseq2_results$Pos)
medz_all = fread(paste0(gtex_dir, 'mpra/all_zscores.txt'))
mpra_zscores = inner_join(deseq2_results, medz_all, by=c('Gene', 'Chr', 'Pos'))
ggplot(mpra_zscores, aes(x=MedZ, y=abs(log2FoldChange_allele))) + geom_point() + theme_bw()

medz_ws_variants = fread('/Users/nicoleferraro/durga_local/data/goats_data/v8_data/mpra/medz_outlier_mpra_watershed_variants.txt') %>%
  group_by(VID,Gene) %>% top_n(1,total_expression_watershed_posterior) %>% sample_n(1) %>% ungroup()

colnames(medz_ws_variants)[11] = 'NEW_VID'
colnames(medz_ws_variants)[7] = 'VID'
deseq2_ws = inner_join(deseq2_results, medz_ws_variants, by=c('Gene','VID'))
deseq2_ws$IsSignificant = sapply(deseq2_ws$padj_allele, function(x) ifelse(x < 0.05, 1, 0))
deseq2_ws$IsSignificant = factor(deseq2_ws$IsSignificant, levels=c(0,1))

ggplot(deseq2_ws, aes(x=total_expression_watershed_posterior, y=abs(log2FoldChange_allele))) +
  geom_point() + theme_bw() +
  annotate("text", x=0.5, y=1.5, label='spearman rho = 0.19, pval = 0.019', cex=5) +
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=16)) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggplot(deseq2_ws, aes(x=total_expression_watershed_posterior, y=-log10(padj_allele))) +
  geom_point() + theme_bw() +
  annotate("text", x=0.5, y=1.5, label='spearman rho = 0.15, pval = 0.063', cex=5) +
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=16)) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

cor.test(deseq2_ws$total_expression_watershed_posterior, -log10(deseq2_ws$padj_allele),method='spearman')

## shendure data
sdata = fread(paste0(gtex_dir, 'mpra/GRCh38_ALL_shendure.tsv'))
colnames(sdata)[1:2] = c('Chr', 'Pos')
sdata$Chr = paste0('chr', sdata$Chr)
medz_outliers = fread(paste0(gtex_dir, 'gtexV8.outliers.v8ciseQTLs.removed.medz.Z2.withrarevars.txt'))
medz_outliers$Pos = as.numeric(medz_outliers$Pos)
medz_outliers = filter(medz_outliers, Chr %in% sdata$Chr)
medz_outliers$Gene = sapply(medz_outliers$Gene, function(x) strsplit(x, '[.]')[[1]][1])
medz_shendure = inner_join(medz_outliers, sdata, by=c('Chr', 'Pos'))



