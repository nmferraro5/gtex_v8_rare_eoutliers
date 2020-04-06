library(data.table)
library(dplyr)
library(ggplot2)
library(bit64)

data_dir = '/Users/nicoleferraro/durga_local/data/goats_data/v8_data/'

ukbb_genomewide_stats = fread(paste0(data_dir, 'coloc/gwas/ukbb.gtex.gwas.all.genomewide.sig.variants.chrpos.txt'))

gtex_gwas = fread(paste0(data_dir, 'coloc/gwas/gtex.traits.beta.gwas.overlap.txt'))
gtex_gwas$Chr = sapply(gtex_gwas$variant, function(x) strsplit(x, ':')[[1]][1])
gtex_gwas$Pos = sapply(gtex_gwas$variant, function(x) strsplit(x, ':')[[1]][2])
gtex_gwas$VID = paste(paste0('chr',gtex_gwas$Chr), gtex_gwas$Pos, sep=':')

trait_table = fread(paste0(data_dir, 'coloc/gwas/ukbb.trait.table.txt'))
colnames(trait_table)[1] = 'Trait'
gtex_gwas = inner_join(gtex_gwas, trait_table, by='Trait')
gtex_gwas = gtex_gwas %>% group_by(Trait) %>% mutate(VPercentile = ntile(abs(beta),100)) %>% ungroup()
gtex_gwas$Pos = as.numeric(gtex_gwas$Pos)
gtex_gwas$Chr = paste0('chr', gtex_gwas$Chr)

get_filtered_data <- function(gtrait) {
  print(gtrait)
  rare_data = filter(gtex_gwas, Trait == gtrait, low_confidence_variant == 'FALSE')
  sig_data = filter(ukbb_genomewide_stats, Trait == gtrait)
  rare_filtered_data_100 = inner_join(rare_data, sig_data, by=c('Chr','Trait')) %>%
    filter(Pos.x > (Pos.y-100000)) %>%
    filter(Pos.x < (Pos.y + 100000)) %>% 
    select(variant,minor_AF,beta,pval,Trait,Chr,Pos.x,VID,Description,VPercentile) %>%
    distinct() %>% mutate(Window = 100000)
  rare_filtered_data_10 = inner_join(rare_data, sig_data, by=c('Chr','Trait')) %>%
    filter(Pos.x > (Pos.y-10000)) %>%
    filter(Pos.x < (Pos.y + 10000)) %>% 
    select(variant,minor_AF,beta,pval,Trait,Chr,Pos.x,VID,Description,VPercentile) %>%
    distinct() %>% mutate(Window = 10000)
  rare_filtered_data_1 = inner_join(rare_data, sig_data, by=c('Chr','Trait')) %>%
    filter(Pos.x > (Pos.y-1000)) %>%
    filter(Pos.x < (Pos.y + 1000)) %>% 
    select(variant,minor_AF,beta,pval,Trait,Chr,Pos.x,VID,Description,VPercentile) %>%
    distinct() %>% mutate(Window = 1000)
  rare_filtered_data = rbind(rare_filtered_data_100, rare_filtered_data_10, rare_filtered_data_1)
  return(rare_filtered_data)
}

gtex_gwas_filtered = do.call(rbind, lapply(unique(gtex_gwas$Trait), get_filtered_data))

te_variants = fread(paste0(data_dir, 'gtexV8.outliers.v8ciseQTLs.globalOutliers.removed.medz.withraresnps.txt'))
ase_variants = fread(paste0(data_dir, 'ASE/gtexV8.outliers.ase.withraresnps.txt'))
sp_variants = fread(paste0(data_dir, 'splicing/gtexV8.outliers.splicing.withraresnps.txt'))
all_variants = unique(c(te_variants$VID, ase_variants$VID, sp_variants$VID))

gtex_gwas_filtered = gtex_gwas_filtered %>% 
  mutate(RareAny = ifelse(VID %in% all_variants, 1, 0))

gtex_gwas_filtered_top = gtex_gwas_filtered %>% group_by(VID,Window) %>% top_n(1,VPercentile) %>% ungroup()

gtex_gwas_filtered_top$RareAny = factor(gtex_gwas_filtered_top$RareAny, levels=c(0,1))
gtex_gwas_filtered_top$Window = factor(gtex_gwas_filtered_top$Window, levels=c(1000, 10000, 100000))
ggplot(gtex_gwas_filtered_top, aes(x=Window,y=VPercentile)) +
  geom_boxplot(width=0.5,aes(fill=RareAny)) +
  labs(fill = "Outlier variant") +
  theme_bw() + xlab('Window around genome-wide significant hit (bp)') + 
  ylab('Best variant effect size percentile') +
  scale_fill_manual(values=c('lightgrey', 'darkcyan')) +
  theme(axis.title = element_text(size=18),
        axis.text = element_text(size=16),
        legend.text=element_text(size=16),
        legend.title=element_text(size=18)) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

wilcox.test(filter(gtex_gwas_filtered_top, Window == 1000, RareAny == 1)$VPercentile,
            filter(gtex_gwas_filtered_top, Window == 1000, RareAny == 0)$VPercentile,
            alternative='g')


eqtl_coloc_results = fread(paste0(data_dir, 'coloc/coloc_sig_results_0.5.txt'))
eqtl_coloc_results$Trait = sapply(eqtl_coloc_results$filename, function(x) strsplit(x, '__PM__')[[1]][1])

eqtl_results = fread(paste0(data_dir, 'coloc/enloc_sig_results_0.5.txt'))
eqtl_results$Trait = sapply(eqtl_results$filename, function(x) strsplit(x, '_w_')[[1]][2])
eqtl_results$Trait = sapply(eqtl_results$Trait, function(x) strsplit(x, '_enloc_output.txt')[[1]][1])
sqtl_results = fread(paste0(data_dir, 'splicing/sqtl_enloc_rcp0.5_collapsed_bygene.txt'))
colnames(sqtl_results) = c('Chr', 'Start', 'End', 'rcp', 'Trait', 'Tissue', 'GeneChr', 'GeneStart', 'GeneEnd', 'Gene')
sqtl_results = sqtl_results %>% select(rcp, Trait, Gene)
