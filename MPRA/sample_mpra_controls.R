library(data.table)
library(dplyr)

data_dir = '/Users/nicoleferraro/durga_local/data/goats_data/v8_data/mpra/mpra_test/'

stop_vars = fread(paste0(data_dir, 'crispr_gtex_stop_candidates.txt'))
stop_exp = melt(fread(paste0(data_dir, 'crispr_stop_candidate_gene_expression.txt')), id.vars=c('Tissue', 'Gene'))

stop_exp_sum = stop_exp %>% group_by(Gene,variable) %>% mutate(NTissues = length(which(abs(value) > 3))) %>% 
  sample_n(1) %>% ungroup()
colnames(stop_exp_sum)[3] = 'Ind'

stop_vars = inner_join(stop_vars, stop_exp_sum, by=c('Gene', 'Ind'))

coloc_genes = fread(paste0(data_dir, 'all_eqtl_coloc_genes.txt'))

stop_vars$IsColoc = sapply(stop_vars$Gene, function(x) ifelse(x %in% coloc_genes$x, 1, 0))

mpra_vars = fread(paste0(data_dir, 'ws5/mpra_gtex_test_candidates_ws5.txt'))
mpra_exp = melt(fread(paste0(data_dir, 'ws5/mpra_candidate_gene_expression.txt')), id.vars=c('Tissue','Gene'))

broad_names = c('Adipose', 'Adipose', 'Adrenal_Gland', rep('Artery',3), rep('Brain', 13), 'Breast', 'Cells_cultured_fibroblasts', 'Cells_EBV', rep('Colon',2), rep('Esophagus',3), rep('Heart',2), 'Kidney', 'Liver', 'Lung', 'Minor_salivary_gland', 'Muscle_Skeletal', 'Nerve_Tibial', 'Ovary', 'Pancreas', 'Pituitary', 'Prostate', rep('Skin', 2), 'Small_Intestine', 'Spleen', 'Stomach', 'Testis', 'Thyroid', 'Uterus', 'Vagina', 'Whole_Blood')
tissue_map = as.data.frame(cbind(unique(mpra_exp$Tissue), broad_names))
colnames(tissue_map)[1] = 'Tissue'
mpra_exp = inner_join(mpra_exp, tissue_map, by='Tissue')

outliers = filter(mpra_exp, abs(value) > 3) %>% group_by(Gene,variable) %>%
  mutate(NTissues = length(unique(broad_names))) %>% sample_n(1) %>% ungroup()

colnames(outliers)[3] = 'Ind'
mpra_vars = merge(mpra_vars, outliers %>% select(-Tissue,-value,-broad_names), by=c('Gene','Ind'),all.x=T)

exp_mpra_vars = filter(mpra_vars, total_expression_watershed_posterior > 0.5, abs(total_expression_outlier_pvalue) < 0.0027)
ase_mpra_vars = filter(mpra_vars, ase_watershed_posterior > 0.5, abs(ase_outlier_pvalue) < 0.0027)

ase_per_tissue = fread(paste0(data_dir, 'ase_per_tissue.txt'))
ase_mpra_vars$GeneID = sapply(ase_mpra_vars$Gene, function(x) strsplit(x, '[.]')[[1]][1])
colnames(ase_per_tissue)[1] = 'GeneID'
colnames(ase_per_tissue)[2] = 'Ind'
ase_mpra_vars = merge(ase_mpra_vars, ase_per_tissue, by=c('GeneID', 'Ind'),all.x=T)

exp_mpra_vars = filter(exp_mpra_vars, NTissues >= 3)
ase_mpra_vars = filter(ase_mpra_vars, NTissues.y >= 3)

ase_mpra_vars = ase_mpra_vars %>% select(-NTissues.x, -GeneID)
colnames(ase_mpra_vars)[17] = 'NTissues'
test_candidates = rbind(exp_mpra_vars, ase_mpra_vars %>% select(-broad_names)) %>% distinct()
test_unique_candidates = test_candidates %>% group_by(Gene,Chr,Pos,ref,allele) %>% top_n(1,MaxWS) %>% sample_n(1) %>% ungroup()

controls = fread(paste0(data_dir, 'ws5/mpra_gtex_test_candidates_controls_ws5_broad.txt'))
controls = controls %>% select(Gene,Chr,Pos,ref,allele,variant_cat,af_gnomad,distTSS) %>% distinct()
controls$af_gnomad = sapply(controls$af_gnomad, function(x) ifelse(is.na(x), 0, x))
test_unique_candidates$af_gnomad = sapply(test_unique_candidates$af_gnomad, function(x) ifelse(is.na(x), 0, x))

get_match <- function(i) {
  tdata = test_unique_candidates[i,]
  tcontrols = filter(controls, Gene == tdata$Gene[1])
  if (nrow(tcontrols) == 0) {
    return(rep(NA, ncol(tcontrols)))
  }
  if (nrow(tcontrols) == 1) {
    tneighbor = get.knnx(tcontrols %>% select(af_gnomad), tdata %>% select(af_gnomad), k=1)
  } else {
    tneighbor = get.knnx(tcontrols %>% select(af_gnomad), tdata %>% select(af_gnomad), k=2)
  }
  tcontrol = tcontrols[tneighbor$nn.index,]
  return(tcontrol)
}

chosen_controls = do.call(rbind, lapply(1:nrow(test_unique_candidates), function(x) get_match(x))) %>%
  filter(!is.na(Gene))

set.seed(kseed) # 9321 was close
sampled_controls = chosen_controls %>% distinct() %>% sample_n(1152)
wilcox.test(test_unique_candidates$af_gnomad, sampled_controls$af_gnomad)
wilcox.test(test_unique_candidates$distTSS, sampled_controls$distTSS)
write.table(sampled_controls, file=paste0(data_dir, 'ws5/mpra_gtex_sampled_controls.txt'), sep='\t', quote=F, row.names=F)

crispr_candidates = fread(paste0(data_dir, 'top_coding_splice_crispr_candidates_ordered.txt'))
crispr_controls = fread(paste0(data_dir, 'crispr_larger_control_pool.txt'))
colnames(crispr_controls)[1:3] = c('Gene', 'Chr', 'Pos')
crispr_candidates$af_gnomad = sapply(as.numeric(crispr_candidates$af_gnomad), function(x) ifelse(is.na(x), 0, x))
crispr_controls$af_gnomad = sapply(as.numeric(crispr_controls$af_gnomad), function(x) ifelse(is.na(x), 0, x))
crispr_candidates$VID = paste(crispr_candidates$Chr, crispr_candidates$Pos, sep=':')
crispr_controls$Chr = paste0('chr', crispr_controls$Chr)
crispr_controls$VID = paste(crispr_controls$Chr, crispr_controls$Pos, sep=':')
crispr_controls = filter(crispr_controls, !(VID %in% crispr_candidates$VID))

get_match <- function(i) {
  tdata = crispr_candidates[i,]
  tcontrols = filter(crispr_controls, Gene == tdata$Gene[1])
  if (nrow(tcontrols) == 0) {
    return(rep(NA, ncol(tcontrols)))
  }
  tneighbor = get.knnx(tcontrols %>% select(Pos,distTSS,af_gnomad), tdata %>% select(Pos,distTSS,af_gnomad), k=1)
  tcontrol = tcontrols[tneighbor$nn.index,]
  return(tcontrol)
}

chosen_controls = do.call(rbind, lapply(1:nrow(crispr_candidates), function(x) get_match(x))) %>%
  filter(!is.na(Gene))

