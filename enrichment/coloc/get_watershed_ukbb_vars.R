library(data.table)
library(dplyr)

data_dir = '/srv/scratch/restricted/GOATs/data_v8/watershed/'
gtex_dir = '/users/nferraro/data/goats_data/v8_data/'

gtex_gwas = fread(paste0(gtex_dir, 'coloc/gwas/gtex.traits.beta.gwas.overlap.txt'))
gtex_gwas$Chr = sapply(gtex_gwas$variant, function(x) strsplit(x, ':')[[1]][1])
gtex_gwas$Pos = sapply(gtex_gwas$variant, function(x) strsplit(x, ':')[[1]][2])
gtex_gwas$VID = paste(paste0('chr',gtex_gwas$Chr), gtex_gwas$Pos, sep=':')

wexp = fread(paste0(data_dir, 'expression.gam.variants.txt'))
wexp$Variant = sapply(wexp$sample_names, function(x) strsplit(x, ':')[[1]][3])
wexp$VID = sapply(wexp$Variant, function(x) paste(strsplit(x, '_')[[1]][1:2],collapse=':'))
wexp = wexp %>% select(-Variant) %>% filter(VID %in% gtex_gwas$VID)

write.table(wexp, file=paste0(data_dir, 'expression.gam.ukbb.variants.txt'), sep='\t',quote=F,row.names=F)

rm(wexp)

wsplice = fread(paste0(data_dir, 'splicing.gam.variants.txt'))
wsplice$Variant = sapply(wsplice$sample_names, function(x) strsplit(x, ':')[[1]][3])
wsplice$VID = sapply(wsplice$Variant, function(x) paste(strsplit(x, '_')[[1]][1:2],collapse=':'))
wsplice = wsplice %>% select(-Variant) %>% filter(VID %in% gtex_gwas$VID)

write.table(wsplice, file=paste0(data_dir, 'splicing.gam.ukbb.variants.txt'), sep='\t',quote=F,row.names=F)

rm(wsplice)

wase = fread(paste0(data_dir, 'ase.gam.variants.txt'))
wase$Variant = sapply(wase$sample_names, function(x) strsplit(x, ':')[[1]][3])
wase$VID = sapply(wase$Variant, function(x) paste(strsplit(x, '_')[[1]][1:2],collapse=':'))
wase = wase %>% select(-Variant) %>% filter(VID %in% gtex_gwas$VID)

write.table(wase, file=paste0(data_dir, 'ase.gam.ukbb.variants.txt'), sep='\t',quote=F,row.names=F)

