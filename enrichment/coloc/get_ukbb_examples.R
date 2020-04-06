library(data.table)
library(dplyr)
library(ggplot2)
library(epitools)
library(scales)
library(qqman)

gtex_gwas = fread(paste0(data_dir, 'gtex.traits.beta.gwas.overlap.txt'))
gtex_gwas$Chr = sapply(gtex_gwas$variant, function(x) strsplit(x, ':')[[1]][1])
gtex_gwas$Pos = sapply(gtex_gwas$variant, function(x) strsplit(x, ':')[[1]][2])
gtex_gwas$VID = paste(paste0('chr',gtex_gwas$Chr), gtex_gwas$Pos, sep=':')
trait_table = fread(paste0(data_dir, 'ukbb.trait.table.txt'))
colnames(trait_table)[1] = 'Trait'
gtex_gwas = inner_join(gtex_gwas, trait_table, by='Trait')
gtex_gwas = gtex_gwas %>% group_by(Trait) %>% mutate(VRank = rank(abs(beta))) %>% ungroup()

exp_ws = fread(paste0(data_dir, 'expression.gam.watershed.ukbb.variants.txt'))
sp_ws = fread(paste0(data_dir, 'splicing.gam.watershed.ukbb.variants.txt'))
ase_ws = fread(paste0(data_dir, 'ase.gam.watershed.ukbb.variants.txt'))

exp_ws = inner_join(exp_ws, gtex_gwas %>% select(minor_AF,low_confidence_variant,beta,pval,Trait,Description,VID, VRank), by='VID')
sp_ws = inner_join(sp_ws, gtex_gwas %>% select(minor_AF,low_confidence_variant,beta,pval,Trait,Description,VID, VRank), by='VID')
ase_ws = inner_join(ase_ws, gtex_gwas %>% select(minor_AF,low_confidence_variant,beta,pval,Trait,Description,VID, VRank), by='VID')

ws_all = rbind(exp_ws %>% mutate(outlier_pvalue = total_expression_outlier_pvalue, gam_posterior = total_expression_gam_posterior, ws_posterior = total_expression_watershed_posterior) %>% 
                 select(-total_expression_outlier_pvalue, -total_expression_gam_posterior, -total_expression_watershed_posterior) %>% mutate(Type = 'TE'),
               sp_ws %>% mutate(outlier_pvalue = splicing_outlier_pvalue, gam_posterior = splicing_gam_posterior, ws_posterior = splicing_watershed_posterior) %>% 
                 select(-splicing_outlier_pvalue, -splicing_gam_posterior, -splicing_watershed_posterior) %>% mutate(Type = 'Splicing'),
               ase_ws %>% mutate(outlier_pvalue = ase_outlier_pvalue, gam_posterior = ase_gam_posterior, ws_posterior = ase_watershed_posterior) %>% 
                 select(-ase_outlier_pvalue, -ase_gam_posterior, -ase_watershed_posterior) %>% mutate(Type = 'ASE'))

rm(exp_ws, sp_ws, ase_ws)

ws_all$Gene = sapply(ws_all$sample_names, function(x) strsplit(x, ':')[[1]][2])
ws_all = filter(ws_all, !is.na(beta))
# choosing 34614 as 75% quantile of ranks
ws_top = filter(ws_all, ws_posterior > 0.9, VRank > 34614, abs(outlier_pvalue) < 0.05, low_confidence_variant == 'FALSE')

## plot examples
cgene = 'ENSG00000177697.17'
gstart = 4709559
gend = 4742043

all_trait_data = fread(paste0(data_dir, 'traits/height.chr9.ukbb.sum.stats.tsv'))
all_trait_data$Chr = sapply(all_trait_data$variant, function(x) strsplit(x, ':')[[1]][1])
all_trait_data$Pos = sapply(all_trait_data$variant, function(x) strsplit(x, ':')[[1]][2])
all_trait_data$VID = paste(paste0('chr', all_trait_data$Chr), all_trait_data$Pos, sep=':')
mp_data = all_trait_data %>% filter(Chr == 9) %>% select(Chr, Pos, pval, variant) %>% filter(!is.na(pval), !is.na(Chr))
mp_data$Pos = as.numeric(mp_data$Pos)
mp_data$Chr = as.numeric(mp_data$Chr)
## manhattan plot
manhattan(mp_data, chr = "Chr", bp = "Pos", p = "pval", snp = "variant", highlight = c('9:4751450:C:T'))

all_trait_plot_data = filter(all_trait_data, Chr == 9, Pos > gstart - 10000, Pos < gend + 10000)
all_trait_plot_data$VOI = sapply(all_trait_plot_data$VID, function(x) ifelse(x == 'chr9:4751450', 'Outlier variant',
                                                                   ifelse(x == 'chr9:4740251', 'lead SNP', 'Other')))
all_trait_plot_data$VOI = factor(all_trait_plot_data$VOI, levels=c('Outlier variant','lead SNP','Other'))
ggplot(all_trait_plot_data %>% filter(low_confidence_variant == 'FALSE'), aes(x=minor_AF,y=abs(beta))) +
  geom_point(aes(color=VOI,size=VOI)) + theme_bw() + guides(size=F) +
  geom_point(data = subset(all_trait_plot_data, VOI != 'Other'),
             aes(color = VOI, size = VOI)) +
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

