library(data.table)
library(dplyr)
library(ggplot2)
require(doMC)
require(foreach)
library(epitools)
library(scales)

registerDoMC(cores = 5)

data_dir = '/users/nferraro/data/goats_data/v8_data/coloc/gwas/'

gtex_gwas = fread(paste0(data_dir, 'gtex.traits.beta.gwas.overlap.txt'))
gtex_gwas$Chr = sapply(gtex_gwas$variant, function(x) strsplit(x, ':')[[1]][1])
gtex_gwas$Pos = sapply(gtex_gwas$variant, function(x) strsplit(x, ':')[[1]][2])
gtex_gwas$VID = paste(paste0('chr',gtex_gwas$Chr), gtex_gwas$Pos, sep=':')
trait_table = fread(paste0(data_dir, 'ukbb.trait.table.txt'))
colnames(trait_table)[1] = 'Trait'
gtex_gwas = inner_join(gtex_gwas, trait_table, by='Trait')

exp_ws = fread(paste0(data_dir, 'expression.watershed.ukbb.variants.txt'))
sp_ws = fread(paste0(data_dir, 'splicing.watershed.ukbb.variants.txt'))
ase_ws = fread(paste0(data_dir, 'ase.watershed.ukbb.variants.txt'))

exp_ws = inner_join(exp_ws, gtex_gwas %>% select(minor_AF,low_confidence_variant,beta,pval,Trait,Description,VID), by='VID')
sp_ws = inner_join(sp_ws, gtex_gwas %>% select(minor_AF,low_confidence_variant,beta,pval,Trait,Description,VID), by='VID')
ase_ws = inner_join(ase_ws, gtex_gwas %>% select(minor_AF,low_confidence_variant,beta,pval,Trait,Description,VID), by='VID')

rm(gtex_gwas)

ws_all = rbind(exp_ws %>% mutate(posterior = total_expression_watershed_posterior) %>% select(-total_expression_watershed_posterior) %>% mutate(Type = 'TE'),
               sp_ws %>% mutate(posterior = splicing_watershed_posterior) %>% select(-splicing_watershed_posterior) %>% mutate(Type = 'Splicing'),
               ase_ws %>% mutate(posterior = ase_watershed_posterior) %>% select(-ase_watershed_posterior) %>% mutate(Type = 'ASE'))

ws_all = ws_all %>% group_by(VID, Type) %>% top_n(1,posterior) %>% ungroup()
ws_all = ws_all %>% group_by(Trait,Type) %>% mutate(VRank = rank(abs(beta))) %>% ungroup()
ws_top = filter(ws_all, posterior > 0.9)
ws_top$Gene = sapply(ws_top$sample_names, function(x) strsplit(x, ':')[[1]][2])

eqtl_coloc = fread('/Users/nicoleferraro/durga_local/data/goats_data/v8_data/coloc/enloc_sig_results_0.5.txt')
sqtl_coloc = fread('/Users/nicoleferraro/durga_local/data/goats_data/v8_data/splicing/sqtl_enloc_rcp0.5_collapsed_bygene.txt')
sqtl_coloc = sqtl_coloc %>% select(V4, V5, V6, V10)
colnames(sqtl_coloc) = c('rcp', 'Trait', 'Tissue', 'Gene')
eqtl_coloc = eqtl_coloc %>% select(gene_id, filename, rcp) %>% mutate(Gene = gene_id) %>% select(-gene_id)

ws_eqtl = inner_join(ws_top, eqtl_coloc, by='Gene')
ws_eqtl$Trait = sapply(ws_eqtl$filename, function(x) strsplit(x, '_w_')[[1]][2])
tkeep = unique(ws_eqtl$Trait)[c(1,4,5,6,9,10,11)]
ws_eqtl = filter(ws_eqtl, Trait %in% tkeep)
ws_sqtl = inner_join(ws_top, sqtl_coloc, by='Gene')
skeep = unique(ws_sqtl$Trait.y)[c(7,8,9)]
ws_sqtl = filter(ws_sqtl, Trait.y %in% skeep)

## plot examples
#cgene = 'ENSG00000177700.5' # for BMI, 837356-842545
cgene = 'ENSG00000143727.15'
gstart = 264393
gend = 278283
all_trait_data = fread(paste0(data_dir, 'height.chr2.ukbb.sum.stats.tsv'))
all_trait_data$Chr = sapply(all_trait_data$variant, function(x) strsplit(x, ':')[[1]][1])
all_trait_data$Pos = sapply(all_trait_data$variant, function(x) strsplit(x, ':')[[1]][2])
all_trait_data$VID = paste(paste0('chr', all_trait_data$Chr), all_trait_data$Pos, sep=':')
all_trait_data = filter(all_trait_data, Pos > gstart - 10000, Pos < gend + 10000)
# eqtl for PUM3: chr9:2825014
# eqtl for TMEM18: chr2:678591
# VID for BMI: chr11:840224, lead coloc snp: chr11_840363 
# VID for height: chr2:281269, lead coloc snp: chr2_274672 (eqtl)
all_trait_data$VOI = sapply(all_trait_data$VID, function(x) ifelse(x == 'chr2:281269', 'Outlier variant',
                                                                     ifelse(x == 'chr2:274672', 'eQTL', 'Other')))
all_trait_data$VOI = factor(all_trait_data$VOI, levels=c('Outlier variant','eQTL','Other'))
#all_height_data = all_height_data %>% filter(abs(beta) < 0.1)
ggplot(all_trait_data %>% filter(low_confidence_variant == 'FALSE'), aes(x=minor_AF,y=abs(beta))) +
  geom_point(aes(color=VOI,size=VOI)) + theme_bw() + guides(size=F) +
  geom_point(data = subset(all_trait_data, VOI != 'Other'),
             aes(color = VOI, size = VOI)) +
  ggtitle('POLR2L in BMI, chr11:840224 with Watershed posterior = 0.99') +
  scale_color_manual(values=c('hotpink', 'aquamarine', 'black')) +
  scale_size_manual(values=c(4,4,2)) +
  theme(axis.text=element_text(size=18),
        axis.title=element_text(size=20),
        legend.title=element_blank(),
        legend.text=element_text(size=18),
        title=element_text(size=20),
        legend.position=c(0.85,0.75)) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ex_data = filter(ws_all, GeneID == cgene, Description == 'Standing height') %>%
  group_by(VID) %>% top_n(1,posterior) %>% sample_n(1) %>% ungroup() %>%
  mutate(HighWS = ifelse(posterior > 0.9, 1, 0))
ex_data$HighWS = factor(ex_data$HighWS, levels=c(0,1))
ggplot(ex_data %>% filter(low_confidence_variant == 'FALSE'), aes(x=minor_AF,y=abs(beta))) + geom_point(aes(color=HighWS))




