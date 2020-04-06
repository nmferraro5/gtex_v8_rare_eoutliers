library(data.table)
library(dplyr)

gwas_dir = '/mnt/lab_data/montgomery/shared/gwas/'
gwas_files = dir(gwas_dir, '*')
eqtl_dir = '/mnt/lab_data/montgomery/shared/datasets/gtex/GTEx_Analysis_2016-01-15_v7/eqtls/all_associations/'
eqtl_files = dir(eqtl_dir, '*.gz')

data_dir = '/users/nferraro/data/goats_data/'
tissues_of_interest = fread(paste0(data_dir, 'coloc/gtex_coloc_tissues_gwas.txt'))$Tissue
eqtl_file_names = sapply(eqtl_files, function(x)
  strsplit(x,'.allpairs.txt.gz')[[1]][1])
kinds = which(eqtl_file_names %in% tissues_of_interest)
eqtl_files_keep = eqtl_files[kinds]

# gwas_of_interest = fread(paste0(data_dir, 'coloc/all_gwas_traits_ranked.txt'), header=F)$V1
# gwas_files_keep = c()
# for (gf in gwas_files) {
#   gname = strsplit(gf, '_')[[1]][1]
#   if (any(grepl(gname, gwas_of_interest))) {
#     gwas_files_keep = c(gwas_files_keep, gf)
#   }
# }
# 
# write.table(gwas_files_keep, file=paste0(data_dir, 'coloc/gwas_dirs_topmed.txt'), sep='\t', quote=F, row.names=F)

# medz_outliers = fread(paste0(data_dir, 'gtexV7.outlier.controls.v7ciseQTLs.removed.medz.txt'), data.table=F)
# knn_outliers = fread(paste0(data_dir, 'gtexV7.outlier.controls.v7ciseQTLs.removed.cor.txt'), data.table=F)
# medz_only_outliers = medz_outliers %>% filter(Y == 'outlier', abs(MedZ) > 2) %>%
#   mutate(Method = 'MEDZ') %>% mutate(Metric = MedZ) %>%
#   select(Ind,Gene,Method,Metric,Y)
# knn_only_outliers = knn_outliers %>% filter(Y == 'outlier') %>%
#   top_n(nrow(medz_only_outliers), -FDR) %>%
#   mutate(Method = 'COR') %>% mutate(Metric = FDR) %>%
#   select(Ind,Gene,Method,Metric,Y)
# 
# rcontrols = filter(medz_outliers, !(Gene %in% medz_only_outliers$Gene)) %>%
#   filter(!(Gene %in% knn_only_outliers$Gene)) %>% distinct(Gene, .keep_all=T) %>%
#   sample_n(4888) %>% mutate(Method = 'Control') %>% mutate(Metric = MedZ) %>%
#   select(Ind,Gene,Method,Metric,Y)
# 
# write.table(rbind(medz_only_outliers, knn_only_outliers, rcontrols), file=paste0(data_dir, 'MEDZ_and_COR_outliers.txt'), sep='\t', quote=F, row.names=F)

medz_cor_outliers = fread(paste0(data_dir, 'MEDZ_and_COR_outliers.txt'))
out_genes = unique(medz_cor_outliers$Gene)

# Only analyze a trait-tissue-gene combo if the gene had at least one eQTL for that 
# gene in that tissue with p-value < 1e-5 and if there was at least one GWAS variant 
# at that locus with p-value < 1e-5.

# do for each eqtl_files_keep
# get_eqtl_data = function(ef) {
#   tissue = strsplit(ef,'.allpairs.txt.gz')[[1]][1]
#   eqtl_data = fread(paste0('zcat ', eqtl_dir, ef)) %>% filter(pval_nominal < 1e-5) %>%
#     filter(gene_id %in% out_genes) %>% mutate(Tissue = tissue)
#   return(eqtl_data)
# }
# 
# all_eqtl_data = do.call(rbind, lapply(eqtl_files_keep, get_eqtl_data))
# out_dir = '/srv/scratch/restricted/GOATs/data_v7/eqtls/'
# write.table(all_eqtl_data, file=paste0(out_dir, 'all_coloc_eqtl_data.txt'),
#             sep='\t', row.names=F, quote=F)
# 
get_gwas_data = function(gdir) {
  gfiles = dir(paste0(gwas_dir, gdir), include.dirs=F, recursive=T)
  print(paste0('NEXT DIRECTORY ', gdir))
  rinds = sapply(gfiles, function(x) ifelse(grepl('zip', x), NA, x))
  gfiles = gfiles[which(gfiles %in% rinds)]
  all_gwases = list()
  for (gf in gfiles) {
    gwas_data = tryCatch({
      if (grepl('gz', gf)) {
        gdata = fread(paste0('zcat ', gwas_dir, gdir, '/', gf), data.table=F, fill=T)
      } else {
        gdata = fread(paste0(gwas_dir, gdir, '/', gf), data.table=F, fill=T)
      }
      pind = which(colnames(gdata) %in% keep_phrases)
      kinds = which(gdata[,pind] < 1e-05)
      gdata = gdata[kinds,]
      gdata$gfile = rep(gf, nrow(gdata))
      gdata$gdir = rep(gdir, nrow(gdata))
      all_gwases[[length(all_gwases)+1]] = gdata
      gdata
    }, error = function(e) {
      NA
    })
  }
  return(all_gwases)
}

keep_phrases = c('pval', 'pvalue', 'P value', 'PVALUE', 'Pvalue', 'P.value', 
                 'P', 'P_top_1_percent', 'p-value', 'p', 'PVAL', 'P-value', 
                 'P_BOLT_LMM_INF', 'pval.GC.SBP','p-value_gc', 'P_VALUE',
                 'CLR_C_BMI_pv','P_FIXED')
gwas_files_keep = fread(paste0(data_dir, 'coloc/gwas_dirs_topmed.txt'))$x
all_gwas_data = lapply(gwas_files_keep, get_gwas_data)
names(all_gwas_data) = gwas_files_keep
saveRDS(all_gwas_data, file='/users/nferraro/data/goats_data/coloc/all_gwas_data.rds')
all_gwas_data = readRDS('/Users/nicoleferraro/durga_local/data/goats_data/coloc/all_gwas_data.rds')
collapsed_gwas = c('SNP', 'CHROM', 'POS', 'PVAL', 'GFILE', 'GDIR')
for (i in seq(1:length(all_gwas_data))) {
  print(names(all_gwas_data)[i])
  print(i)
  gd = all_gwas_data[[i]]
  for (cgd in gd) {
    chrom_ind = which(grepl('chr', colnames(cgd), ignore.case=T),
                      grepl('SNP_hg19', colnames(cgd), ignore.case=T))
    pos_ind = which(grepl('pos', colnames(cgd), ignore.case=T))
    snp_ind = which(grepl('rs', colnames(cgd), ignore.case=T) |
                      grepl('snp', colnames(cgd), ignore.case=T) |
                      grepl('Marker', colnames(cgd), ignore.case=T) |
                      grepl('Chr:Pos', colnames(cgd), ignore.case=T) |
                      grepl('VARIANT', colnames(cgd), ignore.case=T))
    pind = which(colnames(cgd) %in% keep_phrases)
    ninds = which(colnames(cgd) %in% c('gfile', 'gdir'))
    if (length(chrom_ind) == 0 | chrom_ind[1] == snp_ind[1]) { 
      chrom_ind = ncol(cgd) + 1
      cgd$CHROM = rep(NA,nrow(cgd))
    }
    if (length(pos_ind) == 0 | pos_ind[1] == snp_ind[1] | pos_ind[1] == chrom_ind[1]) { 
      pos_ind = ncol(cgd) + 1 
      cgd$POS = rep(NA,nrow(cgd))
    }
    colnames(cgd)[c(snp_ind[1], chrom_ind[1], pos_ind[1], pind, ninds)] = c('SNP', 'CHROM', 'POS', 'PVAL', 'GFILE', 'GDIR')
    collapsed_gwas = rbind(collapsed_gwas, cgd[,c(snp_ind[1], chrom_ind[1], pos_ind[1], pind, ninds)])
  }
}

colnames(collapsed_gwas) = collapsed_gwas[1,]
collapsed_gwas = collapsed_gwas[-1,]
write.table(collapsed_gwas, file='/Users/nicoleferraro/durga_local/data/goats_data/coloc/collapsed_gwas_stats.txt', sep='\t', quote=F, row.names=F)

collapsed_gwas = fread('/Users/nicoleferraro/durga_local/data/goats_data/coloc/collapsed_gwas_stats.txt')
noinfo_inds = which(is.na(collapsed_gwas$POS) & (grepl('rs', collapsed_gwas$SNP)))
chinds = which((is.na(collapsed_gwas$CHROM) & !(grepl('rs', collapsed_gwas$SNP))) |
                 (is.na(collapsed_gwas$POS) & !(grepl('rs', collapsed_gwas$SNP))))
chrom_info = sapply(chinds, function(x)
  strsplit(collapsed_gwas$SNP[x], ':')[[1]][1])
pos_info = sapply(chinds, function(x)
  strsplit(collapsed_gwas$SNP[x], ':')[[1]][2])
collapsed_gwas$CHROM[chinds] = unlist(chrom_info)
collapsed_gwas$POS[chinds] = unlist(pos_info)

ochinds = which(grepl(':', collapsed_gwas$CHROM))
chrom_info = sapply(ochinds, function(x)
  strsplit(collapsed_gwas$CHROM[x], ':')[[1]][1])
pos_info = sapply(ochinds, function(x)
  strsplit(collapsed_gwas$CHROM[x], ':')[[1]][2])
collapsed_gwas$CHROM[ochinds] = unlist(chrom_info)
collapsed_gwas$POS[ochinds] = unlist(pos_info)

collapsed_gwas_nopos = collapsed_gwas[noinfo_inds,]
write.table(collapsed_gwas_nopos$SNP, file='/Users/nicoleferraro/durga_local/data/goats_data/coloc/rsid_no_pos.txt', sep='\t', quote=F, row.names=F, col.names=F)
# Ran the below with the dbSNP file on durga to get unknown positions - discrepancies though with reported chom/pos from gwas files
# zgrep -wFf /users/nferraro/data/goats_data/coloc/rsid_no_pos.txt hg19_snp150.txt.gz > /users/nferraro/data/goats_data/coloc/rsid_pos.txt
# rsid_pos = fread('/Users/nicoleferraro/durga_local/data/goats_data/coloc/rsid_pos_hg38.txt', data.table=F)[,c(2:5)]
# colnames(rsid_pos) = c('CHROM', 'START', 'STOP', 'RSID')
# gwas_chroms = sapply(collapsed_gwas_nopos$SNP, function(x)
#   filter(rsid_pos, RSID == x)$CHROM)
# gwas_pos = sapply(collapsed_gwas_nopos$SNP, function(x)
#   filter(rsid_pos, RSID == x)$START)

cinds = which(grepl('chr', collapsed_gwas$CHROM, ignore.case=T))
collapsed_gwas$CHROM[cinds] = sapply(collapsed_gwas$CHROM[cinds], function(x)
  substr(x, start=4, stop=nchar(x)))
collapsed_gwas$CHROM = sapply(collapsed_gwas$CHROM, function(x)
  ifelse(x == '05', 5, x))
eqtl_data = fread('/Users/nicoleferraro/durga_local/data/goats_data/coloc/all_coloc_eqtl_data.txt')
echroms = sapply(eqtl_data$variant_id, function(x)
  strsplit(x, '_')[[1]][1])
epos = sapply(eqtl_data$variant_id, function(x)
  strsplit(x, '_')[[1]][2])
eqtl_data$CHROM = echroms
eqtl_data$POS = epos

eqtl_data = eqtl_data %>% mutate(CPOS = paste(CHROM,POS,sep='_'))
collapsed_gwas = collapsed_gwas %>% mutate(CPOS = paste(CHROM,POS,sep='_'))

eqtl_gwas = merge(eqtl_data, collapsed_gwas, by='CPOS') %>%
  select(CHROM.x, POS.x, gene_id, variant_id, pval_nominal,Tissue,PVAL,GFILE,GDIR) %>%
  group_by(CHROM.x,POS.x,Tissue,GFILE) %>% sample_n(1) %>% ungroup()
colnames(eqtl_gwas) = c('CHROM', 'POS', 'GENE', 'VARIANT', 'EQTL_PVAL','TISSUE', 'GWAS_PVAL', 'GWASFILE', 'GWASDIR')
gwas_metadata = fread('/Users/nicoleferraro/durga_local/data/goats_data/coloc/gwas_metadata.txt')

write.table(eqtl_gwas, file='/Users/nicoleferraro/durga_local/data/goats_data/coloc/coloc_candidates.txt',
            sep='\t', quote=F, row.names=F)

