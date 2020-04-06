library(data.table)
library(dplyr)
library(ggplot2)

data_dir = '/Users/nicoleferraro/durga_local/data/goats_data/coloc/'
coloc_data = fread(paste0(data_dir, 'all_colocalization_tests_v2.txt'))
coloc_data = coloc_data[-9833,]
coloc_data$clpp = as.numeric(coloc_data$clpp)
coloc_data$clpp_mod = as.numeric(coloc_data$clpp_mod)

# Calculate significance of overlap
m = #number of white balls in the urn

# Compare to v7 outliers
outlier_dir = '/Users/nicoleferraro/durga_local/data/goats_data/v8_data/'
medz.outliers = fread(paste0(outlier_dir, 'gtexV8.outliers.v8ciseQTLs.removed.medz.txt')) %>%
  filter(abs(MedZ) > 2)
cor.outliers = fread(paste0(outlier_dir, 'gtexV8.outliers.v8ciseQTLs.removed.knn.txt')) %>%
  top_n(nrow(medz.outliers),-FDR)
stfz.outliers = fread(paste0(outlier_dir, 'gtexV8.outliers.v8ciseQTLs.removed.stfz.txt')) %>%
  top_n(nrow(medz.outliers),-FDR)

medz.outlier.genes = sapply(medz.outliers$Gene, function(x) strsplit(x,'[.]')[[1]][1])
cor.outlier.genes = sapply(cor.outliers$Gene, function(x) strsplit(x,'[.]')[[1]][1])
stfz.outlier.genes = sapply(stfz.outliers$Gene, function(x) strsplit(x,'[.]')[[1]][1])
all.genes = unique(c(medz.outlier.genes, cor.outlier.genes, stfz.outlier.genes))
coloc.genes = sapply(coloc_data$feature, function(x) strsplit(x, '[.]')[[1]][1])
m = length(which(all.genes %in% coloc.genes))
n = length(unique(coloc.genes)) - m
k = length(unique(filter(coloc_data, clpp_mod > 0.3)$feature))
q = 453
phyper(q,m,n,k,lower.tail=FALSE)



coloc_data$feature = sapply(coloc_data$feature, function(x) strsplit(x,'[.]')[[1]][1])

coloc_data = coloc_data %>% 
  mutate(MedzOutlier = ifelse(feature %in% medz.outlier.genes, 1, 0)) %>%
  mutate(StfzOutlier = ifelse(feature %in% stfz.outlier.genes, 1, 0)) %>%
  mutate(CorOutlier = ifelse(feature %in% cor.outlier.genes, 1, 0))
coloc_data$MedzOutlier = factor(coloc_data$MedzOutlier, levels=c(0,1))
coloc_data$StfzOutlier = factor(coloc_data$StfzOutlier, levels=c(0,1))
coloc_data$CorOutlier = factor(coloc_data$CorOutlier, levels=c(0,1))

clppplot = ggplot(coloc_data, aes(x=StfzOutlier, y=clpp_mod)) +
  geom_boxplot()
clppplot

coloc_medz = filter(coloc_data, MedzOutlier == 1, clpp_mod > 0.3)
coloc_stfz = filter(coloc_data, StfzOutlier == 1, clpp_mod > 0.3)
coloc_cor = filter(coloc_data, CorOutlier == 1, clpp_mod > 0.3)

all_V7coloc = rbind(coloc_medz %>% mutate(Method='MEDZ'),
                        coloc_stfz %>% mutate(Method='STFZ'),
                        coloc_cor %>% mutate(Method='COR'))

all_V7outliers = cbind(rbind(medz.outliers[,1:2], stfz.outliers[,1:2], cor.outliers[,1:2]),
                       c(rep('MEDZ',4888), rep('STFZ',4888), rep('COR',4888)))
colnames(all_V7outliers) = c('Ind','Gene','Method')

all_V8coloc = rbind(coloc_medz %>% mutate(Method='MEDZ'),
                    coloc_stfz %>% mutate(Method='STFZ'),
                    coloc_cor %>% mutate(Method='COR'))

all_V8outliers = cbind(rbind(medz.outliers[,1:2], stfz.outliers[,1:2], cor.outliers[,1:2]),
                       c(rep('MEDZ',6215), rep('STFZ',6215), rep('COR',6215)))
colnames(all_V8outliers) = c('Ind','Gene','Method')

all_V7outliers$Gene = sapply(all_V7outliers$Gene, function(x) strsplit(x,'[.]')[[1]][1])
all_V8outliers$Gene = sapply(all_V8outliers$Gene, function(x) strsplit(x,'[.]')[[1]][1])

colnames(all_V7coloc)[4] = 'Gene'
all_V7 = merge(all_V7coloc, all_V7outliers, by='Gene')
colnames(all_V8coloc)[4] = 'Gene'
all_V8 = merge(all_V8coloc, all_V8outliers, by='Gene')

write.table(all_V8, file='/Users/nicoleferraro/durga_local/data/goats_data/coloc/coloc_V8outliers.txt', sep='\t', quote=F, row.names=F)

## Old analysis
# data_dir = '/Users/nicoleferraro/durga_local/data/goats_data/coloc/'
# # get list of all possible gwas traits
# traits = fread(paste0(data_dir, 'all_gwas.txt'), header=F)
# traits$V1 = sapply(traits$V1, function(x)
#   strsplit(x, '_')[[1]][1])
# uniq_traits = unique(traits$V1)
# 
# coloc_files = dir(data_dir, '*_status.txt', full=T)
# 
# outliers = fread(paste0(data_dir, 'tissue_specific_cor_outliers.txt')) %>%
#   mutate(GT = paste(Gene,Specific.Group,sep='_'))
# coloc_data = do.call(rbind, lapply(coloc_files, function(x) fread(x,data.table=F))) %>%
#   mutate(Tissue = substr(eqtl_file, start=1, stop=nchar(eqtl_file)-16)) %>%
#   mutate(GT = paste(feature,Tissue,sep='_'))
# 
# sig_colocs = filter(coloc_data, clpp_mod > 0.05)
# sig_ts_colocs = merge(sig_colocs, outliers, by='GT')
# nrow(sig_ts_colocs)
# 
# coloc_outliers = merge(coloc_data, outliers, by='GT')
# noncoloc_outliers = filter(coloc_data, feature %in% coloc_outliers$feature,
#                            Tissue %in% coloc_outliers$Tissue,
#                            !(GT %in% coloc_outliers$GT))
# 
# all_clpps = as.data.frame(cbind(c(coloc_outliers$clpp_mod, noncoloc_outliers$clpp_mod),
#                   c(rep('TS_Outlier', nrow(coloc_outliers)),
#                     rep('NonTS', nrow(noncoloc_outliers)))))
# 
# all_clpps$V1 = as.numeric(as.character(all_clpps$V1))
# all_clpps = all_clpps %>% mutate(Status = ifelse(V2 == 'NonTS', 0, 1))
# cplot = ggplot(all_clpps, aes(x=V2,y=V1)) + geom_boxplot() + theme_bw() +
#   xlab('') + ylab('Modified posterior probability') +
#   theme(axis.text=element_text(size=16),
#         axis.title=element_text(size=18)) +
#   annotate("text", x=1.5, y=0.125, label="p = 0.27", cex=5)
# cplot
# 
# summary(lm(V1~Status, data=all_clpps))
