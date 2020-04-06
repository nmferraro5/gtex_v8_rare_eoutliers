library(data.table)
library(dplyr)
library(ggplot2)

data_dir = '/users/nferraro/data/goats_data/v8_data/'

## Intersect with outliers
medz_outliers = fread(paste0(data_dir, 'gtexV8.outliers.v8ciseQTLs.removed.medz.Z2.txt')) %>%
  filter(abs(MedZ) > 3)

sp_outliers = fread(paste0(data_dir, 'splicing/gtexV8.outliers.v8ciseQTLs.removed.splicing.Z2.txt')) %>%
  filter(value < 0.0027)

ase_outliers = fread(paste0(data_dir, 'ASE/median.ad.scores.uncorrected.no.global.outliers.v8.vg.outliers.only.tsv')) %>%
  filter(value < 0.0027)
colnames(ase_outliers)[1] = 'Gene'
gene_names = fread(paste0(data_dir, 'gencode.v26.GRCh38.genes.bed')) %>%
  select(-V1,-V2,-V3) %>% mutate(Gene = V4)
gene_names$Gene = sapply(gene_names$Gene, function(x) strsplit(x, '[.]')[[1]][1])
ase_outliers = inner_join(ase_outliers, gene_names, by='Gene') %>%
  mutate(GID = paste(variable,V4, sep='_'))

all_outliers = rbind(medz_outliers %>% mutate(value = MedZ) %>% select(-MedZ) %>% mutate(Type = 'TE'),
                     sp_outliers %>%  mutate(Type = 'Splicing'),
                     ase_outliers %>% mutate(Ind = variable,Gene = V4) %>% select(Gene,Ind,value) %>% mutate(Type = 'ASE'))

genes_with_rare = fread(paste0(data_dir, 'all_genes_with_rare_variants_combined.txt'),header=F)
colnames(genes_with_rare) = c('Ind','Gene')
genes_with_rare$HasRare = rep(1,nrow(genes_with_rare))
all_outliers = merge(all_outliers, genes_with_rare, by=c('Ind', 'Gene'), all.x=T)
all_outliers$HasRare = sapply(all_outliers$HasRare, function(x) ifelse(is.na(x), 0, 1))
all_outliers = filter(all_outliers, Ind %in% genes_with_rare$Ind)

watershed_data = fread(paste0(data_dir, 'watershed.variants.0.5.v8.update.txt'))
colnames(watershed_data)[3] = 'Type'
ws_0.5_count = watershed_data %>% group_by(Type, Ind, Gene) %>% top_n(1,posterior) %>% sample_n(1) %>%
  ungroup() %>% group_by(Type,Ind) %>% mutate(NW5 = n()) %>% 
  sample_n(1) %>% ungroup() %>% select(Ind,Gene,NW5,Type)

ws_0.9_count = watershed_data %>% filter(posterior > 0.9) %>%
  group_by(Type, Ind, Gene) %>% top_n(1,posterior) %>% sample_n(1) %>%
  ungroup() %>% group_by(Type,Ind) %>% mutate(NW9 = n()) %>% 
  sample_n(1) %>% ungroup() %>% select(Ind,Gene,NW9,Type)

ws_0.7_count = watershed_data %>% filter(posterior > 0.7) %>%
  group_by(Type, Ind, Gene) %>% top_n(1,posterior) %>% sample_n(1) %>%
  ungroup() %>% group_by(Type,Ind) %>% mutate(NW7 = n()) %>% 
  sample_n(1) %>% ungroup() %>% select(Ind,Gene,NW7,Type)

# all_outliers = merge(all_outliers, watershed_data, by=c('Ind','Gene','Type'),all.x=T)
# 
# all_outliers = all_outliers %>% group_by(Ind,Gene,Type) %>% 
#   mutate(MaxP = max(posterior,na.rm=T)) %>% sample_n(1) %>% ungroup()
# all_outliers$MaxP = sapply(all_outliers$MaxP, function(x) ifelse(is.infinite(x), NA, x))

## integrate coloc data
eqtl_coloc_results = fread(paste0(data_dir, 'coloc/coloc_sig_results_0.5.txt'))
eqtl_coloc_results$Trait = sapply(eqtl_coloc_results$filename, function(x) strsplit(x, '__PM__')[[1]][1])
eqtl_results = fread(paste0(data_dir, 'coloc/enloc_sig_results_0.5.txt'))
eqtl_results$Trait = sapply(eqtl_results$filename, function(x) strsplit(x, '_w_')[[1]][2])
eqtl_results$Trait = sapply(eqtl_results$Trait, function(x) strsplit(x, '_enloc_output.txt')[[1]][1])
sqtl_results = fread(paste0(data_dir, 'splicing/sqtl_enloc_rcp0.5_collapsed_bygene.txt'))
colnames(sqtl_results) = c('Chr', 'Start', 'End', 'rcp', 'Trait', 'Tissue', 'GeneChr', 'GeneStart', 'GeneEnd', 'Gene')
sqtl_results = sqtl_results %>% select(rcp, Trait, Gene)

coloc_genes = unique(c(eqtl_coloc_results$gene, eqtl_results$gene_id, sqtl_results$Gene))
all_outliers$IsRareColoc = sapply(1:nrow(all_outliers), function(x)
  ifelse(all_outliers$HasRare[x] == 1 && all_outliers$Gene[x] %in% coloc_genes, 1, 0))

# sum_outliers = all_outliers %>% group_by(Ind,Type) %>%
#   mutate(NG = length(unique(Gene))) %>%
#   mutate(NRV = sum(HasRare)) %>%
#   mutate(NCRV = sum(IsRareColoc)) %>%
#   mutate(NWV = length(which(!is.na(MaxP)))) %>% sample_n(1) %>% ungroup()

sum_outliers = all_outliers %>% group_by(Ind,Type) %>%
  mutate(NG = length(unique(Gene))) %>%
  mutate(NRV = sum(HasRare)) %>%
  mutate(NCRV = sum(IsRareColoc)) %>% sample_n(1) %>% ungroup()

sum_outliers = merge(sum_outliers, ws_0.5_count %>% select(-Gene), by=c('Ind', 'Type'), all.x=T)
sum_outliers = merge(sum_outliers, ws_0.9_count %>% select(-Gene), by=c('Ind', 'Type'), all.x=T)
sum_outliers = merge(sum_outliers, ws_0.7_count %>% select(-Gene), by=c('Ind', 'Type'), all.x=T)


sum_melted = melt(sum_outliers %>% select(Ind,Type,NG,NRV,NCRV,NW5,NW7,NW9), id.vars=c('Ind','Type'))

sum_melted$variable = sapply(sum_melted$variable, function(x)
  ifelse(x == 'NG', 'Outliers',
         ifelse(x == 'NRV', 'Outliers w/ RV', 
                ifelse(x == 'NCRV', 'Coloc Outliers w/ RV', 
                       ifelse(x == 'NW5', 'Outliers w/ high WS > 0.5', 
                              ifelse(x == 'NW7', 'Outliers w/ high WS > 0.7', 'Outliers w/ high WS > 0.9'))))))
# sum_melted$variable = factor(sum_melted$variable, levels=c('Outliers', 'Outliers w/ RV', 'Coloc Outliers w/ RV', 'Outliers w/ high WS'))
# 
# splot = ggplot(sum_melted, aes(x=variable,y=value)) + geom_boxplot() + theme_bw() +
#   xlab('') + ylab('Number per individual') +
#   theme(axis.text=element_text(size=14),
#         axis.title=element_text(size=18),
#         strip.text.y = element_text(size=18)) +
#   facet_grid(Type~.,scales='free')

save(sum_melted, file=paste0(data_dir, 'gtexV8.outlier.Z3.variant.count.per.individual.v8.update.ws.RData'))








