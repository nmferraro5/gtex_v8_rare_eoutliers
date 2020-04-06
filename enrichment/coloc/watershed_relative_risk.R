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


ws_all_top = ws_all %>% filter(!is.na(outlier_pvalue)) %>% filter(!is.na(beta)) %>% group_by(VID,Trait,Type) %>% top_n(1,ws_posterior) %>%
  sample_n(1) %>% ungroup() %>% group_by(Trait,Type) %>% 
  mutate(pbin = ifelse(ws_posterior > 0.5, 'High', 'Low')) %>% ungroup() %>% 
  mutate(posterior = ws_posterior) %>%
  select(-gam_posterior, -outlier_pvalue, -ws_posterior)
gam_all_top = ws_all %>% filter(!is.na(outlier_pvalue)) %>% filter(!is.na(beta)) %>% group_by(VID,Trait,Type) %>% top_n(1,gam_posterior) %>%
  sample_n(1) %>% ungroup() %>% group_by(Trait,Type) %>% 
  mutate(pbin = ifelse(gam_posterior > 0.025, 'High', 'Low')) %>% ungroup() %>% 
  mutate(posterior = gam_posterior) %>%
  select(-gam_posterior, -outlier_pvalue, -ws_posterior)

all_top = rbind(ws_all_top %>% mutate(Model = 'Watershed'),
                gam_all_top %>% mutate(Model = 'GAM'))
#all_top$Gene = sapply(all_top$sample_names, function(x) strsplit(x, ':')[[1]][2])
all_top$pbin = factor(all_top$pbin, levels=c('Low', 'High'))

#kws_genes = filter(all_top, Model == 'Watershed', pbin == 1)$Gene

## top 100 variants - proportion in each bin
top_prop = all_top %>% group_by(Model,Type,Trait) %>% top_n(258,abs(beta)) %>% ungroup() %>%
  group_by(VID,Type) %>% sample_n(1) %>% ungroup() %>%
  group_by(Model,Type) %>% mutate(NV = n()) %>% ungroup() %>%
  group_by(Model,Type,Trait,pbin) %>% mutate(PropBin = n()/NV) %>% sample_n(1) %>% ungroup() %>%
  select(Model,Type,Trait,pbin,PropBin) %>% mutate(Cat='Top 1%')
background_prop = all_top %>% group_by(Model,Type) %>% mutate(NV = n()) %>% ungroup() %>%
  group_by(VID,Type) %>% sample_n(1) %>% ungroup() %>%
  group_by(Model,Type,Trait,pbin) %>% mutate(PropBin = n()/NV) %>% sample_n(1) %>% ungroup() %>%
  select(Model,Type,Trait,pbin,PropBin) %>% mutate(Cat = 'Background')

top_plot = rbind(top_prop, background_prop) %>%
  group_by(pbin,Model,Type,Cat) %>% sample_n(1) %>% ungroup()
ggplot(top_plot %>% filter(pbin == 'High', Model == 'Watershed'), aes(x=Type,y=PropBin)) + 
  geom_bar(aes(fill=Cat), stat='identity', position='dodge') +
  ylab('Proportion of variants with high posterior') + xlab('') +
  theme_bw() + scale_fill_manual(values=c('grey', 'firebrick4')) +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=17),
        legend.title=element_blank(),
        legend.text=element_text(size=15),
        strip.text.x=element_text(size=15)) 

ws_risks = data.frame(Riskratio = numeric(), Lower = numeric(), Upper = numeric(), Pval = numeric(), Type = character())
for (gtype in unique(all_top$Type)) {
  top_data = filter(all_top, Type == gtype, Model == 'Watershed') %>% 
    group_by(Model,Trait) %>% top_n(461,abs(beta)) %>% ungroup() %>%
    group_by(VID) %>% sample_n(1) %>% ungroup()
  other_data = filter(all_top, Type == gtype, Model == 'Watershed', !(VID %in% top_data$VID)) %>%
    group_by(VID) %>% sample_n(1) %>% ungroup()
  no_no = nrow(filter(other_data, pbin == 'Low'))
  no_yes = nrow(filter(other_data, pbin == 'High'))
  yes_no = nrow(filter(top_data, pbin == 'Low'))
  yes_yes = nrow(filter(top_data, pbin == 'High'))
  counttable = rbind(c(no_no, no_yes), c(yes_no, yes_yes))
  rr = epitab(counttable,method="riskratio")$tab
  ws_risks = rbind(ws_risks, data.frame(Riskratio = rr[2,5], Lower = rr[2,6], Upper = rr[2,7], Pval = rr[2,8], Type = gtype))
}

gam_risks = data.frame(Riskratio = numeric(), Lower = numeric(), Upper = numeric(), Pval = numeric(), Type = character())
for (gtype in unique(all_top$Type)) {
  top_data = filter(all_top, Type == gtype, Model == 'GAM') %>% 
    group_by(Model,Trait) %>% top_n(461,abs(beta)) %>% ungroup() %>%
    group_by(VID) %>% sample_n(1) %>% ungroup()
  other_data = filter(all_top, Type == gtype, Model == 'GAM', !(VID %in% top_data$VID)) %>%
    group_by(VID) %>% sample_n(1) %>% ungroup()
  no_no = nrow(filter(other_data, pbin == 'Low'))
  no_yes = nrow(filter(other_data, pbin == 'High'))
  yes_no = nrow(filter(top_data, pbin == 'Low'))
  yes_yes = nrow(filter(top_data, pbin == 'High'))
  counttable = rbind(c(no_no, no_yes), c(yes_no, yes_yes))
  rr = epitab(counttable,method="riskratio")$tab
  gam_risks = rbind(gam_risks, data.frame(Riskratio = rr[2,5], Lower = rr[2,6], Upper = rr[2,7], Pval = rr[2,8], Type = gtype))
}

both_risks = rbind(ws_risks %>% mutate(Model = 'Watershed'),
                   gam_risks %>% mutate(Model = 'GAM'))
ggplot(both_risks %>% filter(Riskratio < 10), aes(x=Type,y=Riskratio,Group=Model)) + 
  geom_point(aes(color=Model),size=4,position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width=0,position=position_dodge(width=0.5)) +
  geom_hline(yintercept=1,color='grey') + xlab('') + 
  scale_color_manual(values=c('grey','firebrick4')) +
  ylab('Relative risk of high posterior variant in top 1% effect size traits') +
  theme_bw() + theme(axis.text=element_text(size=16),
                     axis.title=element_text(size=18),
                     legend.title=element_blank(),
                     legend.text=element_text(size=16),
                     legend.position=c(0.5,0.85))
  
  
  
top_data = filter(all_top, Model == 'Watershed') %>% 
  group_by(Model,Trait) %>% top_n(2307,abs(beta)) %>% ungroup() %>%
  group_by(VID) %>% top_n(1,posterior) %>% sample_n(1) %>% ungroup()
other_data = filter(all_top, Model == 'Watershed', !(VID %in% top_data$VID)) %>%
  group_by(VID) %>% top_n(1,posterior) %>% sample_n(1) %>% ungroup()
no_no = nrow(filter(other_data, pbin == 'Low'))
no_yes = nrow(filter(other_data, pbin == 'High'))
yes_no = nrow(filter(top_data, pbin == 'Low'))
yes_yes = nrow(filter(top_data, pbin == 'High'))
counttable = rbind(c(no_no, no_yes), c(yes_no, yes_yes))
rr = epitab(counttable,method="riskratio")$tab
ws_risks = rbind(ws_risks, data.frame(Riskratio = rr[2,5], Lower = rr[2,6], Upper = rr[2,7], Pval = rr[2,8], Type = 'All'))

top_data = filter(all_top, Model == 'GAM') %>% 
  group_by(Model,Trait) %>% top_n(2307,abs(beta)) %>% ungroup() %>%
  group_by(VID) %>% top_n(1,posterior) %>% sample_n(1) %>% ungroup()
other_data = filter(all_top, Model == 'GAM', !(VID %in% top_data$VID)) %>%
  group_by(VID) %>% top_n(1,posterior) %>% sample_n(1) %>% ungroup()
no_no = nrow(filter(other_data, pbin == 'Low'))
no_yes = nrow(filter(other_data, pbin == 'High'))
yes_no = nrow(filter(top_data, pbin == 'Low'))
yes_yes = nrow(filter(top_data, pbin == 'High'))
counttable = rbind(c(no_no, no_yes), c(yes_no, yes_yes))
rr = epitab(counttable,method="riskratio")$tab
gam_risks = rbind(gam_risks, data.frame(Riskratio = rr[2,5], Lower = rr[2,6], Upper = rr[2,7], Pval = rr[2,8], Type = 'All'))

both_risks = rbind(ws_risks %>% mutate(Model = 'Watershed'),
                   gam_risks %>% mutate(Model = 'GAM'))
ggplot(both_risks %>% filter(Type == 'All'), aes(x=Model,y=Riskratio)) + 
  geom_point(size=4) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width=0) +
  geom_hline(yintercept=1,color='grey') + xlab('') + 
  scale_color_manual(values=c('grey','firebrick4')) +
  ylab('Relative risk of high posterior variant in top 1% effect size traits') +
  theme_bw() + theme(axis.text=element_text(size=16),
                     axis.title=element_text(size=18))



