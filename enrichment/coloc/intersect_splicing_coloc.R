library(data.table)
library(dplyr)

data_dir = '/Users/nicoleferraro/durga_local/data/goats_data/v8_data/'

sqtl_data = fread(paste0(data_dir, 'splicing/sqtl_enloc_rcp0.5_collapsed_update.txt'))
sqtl_data$Chr = sapply(sqtl_data$molecular_qtl_trait, function(x) paste0('chr', strsplit(x, '_')[[1]][2]))
sqtl_data$Start = sapply(sqtl_data$molecular_qtl_trait, function(x) strsplit(x, '_')[[1]][3])
sqtl_data$End = sapply(sqtl_data$molecular_qtl_trait, function(x) strsplit(x, '_')[[1]][4])

sqtl_data_out = sqtl_data %>% select(Chr,Start,End,locus_rcp,Trait,Tissue)

write.table(sqtl_data_out, file=paste0(data_dir, 'sqtl_enloc_rcp0.5_collapsed_bygene.txt'),sep='\t',quote=F,row.names=F)

sqtl_data = fread(paste0(data_dir, 'splicing/sqtl_enloc_rcp0.5_collapsed_bygene.txt'))
sqtl_data = sqtl_data %>% select(-V7,-V8,-V9)
colnames(sqtl_data) = c('Chr','Start','End','rcp','Trait','Tissue','Gene')