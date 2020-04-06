#!/usr/bin/env Rscript

# Calculate absolute risk of outlier given RV

## Load required packages
require(data.table)
require(dplyr)
library(ggplot2)
require(foreach)
require(doMC)

registerDoMC(cores = 3)

data_dir = '/users/nferraro/data/goats_data/v8_data/'
goats_dir = Sys.getenv('RAREDIR')

print('Reading in outliers')
zthresh = 3
medz_data = fread(paste0(data_dir, 'gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.medz.txt'))
#medz_data$GeneID = sapply(medz_data$Gene, function(x) strsplit(x,'[.]')[[1]][1])
medz_data = medz_data %>% mutate(UID = paste(Ind,Gene,sep='_'))
medz_outliers = medz_data %>% filter(abs(MedZ) > zthresh)
splicing_data = fread(paste0(goats_dir, '/data_v8/splicing/cross_tissue_covariate_method_none_no_global_outliers_ea_only_emperical_pvalue_gene_level.txt'))
splicing_data = melt(splicing_data) 
#splicing_data$GeneID = sapply(splicing_data$CLUSTER_ID, function(x) strsplit(x,'[.]')[[1]][1])
splicing_data = splicing_data %>% mutate(UID = paste(variable,CLUSTER_ID,sep='_'))
splicing_outliers = filter(splicing_data,!is.nan(value),value < 2*pnorm(-abs(zthresh)))
ase_data = melt(fread(paste0(goats_dir, '/data_v8/ase/combined.ad.scores.in.MEDIAN_4_10_update.tsv')))
colnames(ase_data) = c('GeneID', 'SampleName', 'DOT.score')
gene.mapper = fread(paste0(data_dir, 'gencode.v26.GRCh38.genes.bed'), header=F)
gene.mapper$GeneID = sapply(gene.mapper$V4, function(x) strsplit(x, '[.]')[[1]][1])
colnames(gene.mapper)[4] = 'Gene'
ase_data = inner_join(ase_data, gene.mapper %>% select(Gene,GeneID), by='GeneID')
ase_data = ase_data %>% mutate(UID = paste(SampleName,Gene,sep='_'))
ase_outliers = ase_data %>% filter(DOT.score < 2*pnorm(-abs(zthresh)))

print('Reading in variant annotations')
variant_data = fread(paste0('zcat ', goats_dir, '/data_v8/variant_annotations/gtex_v8_rare_GxIxV_variant_annotations.txt.gz'),header=T,fill=T)
colnames(variant_data)[5] = 'tier2'
#variant_data$gene_id = as.character(variant_data$gene_id)
#variant_data$variant_cat = as.character(variant_data$variant_cat)
#variant_data$tier2 = as.character(variant_data$tier2)

new_cats = sapply(1:nrow(variant_data), function(x) ifelse(variant_data$variant_cat[x] == 'splice', variant_data$tier2[x], variant_data$variant_cat[x]))
variant_data$variant_cat = new_cats
all_cats = c('no_variant','other_noncoding','coding','conserved_noncoding','TSS','stop','frameshift','splice_region_variant','splice_donor_variant','splice_acceptor_variant','TE','INV','BND','DEL','CNV','DUP')
variant_data = filter(variant_data, variant_cat %in% all_cats, sv_v7 == 1) %>% select(indiv_id, gene_id, variant_cat)
variant_data$indiv_id = paste0('GTEX-', variant_data$indiv_id)
variant_data$UID = paste(variant_data$indiv_id, variant_data$gene_id, sep='_')

print('Starting absolute risk calculation')
random_abs_risk = foreach(i=1:1000, .combine = rbind) %dopar% {
    print(paste0('Iteration ', i))
    random_medz_uids = sample(medz_data$UID, nrow(medz_outliers))
    random_sp_uids = sample(splicing_data$UID, nrow(splicing_outliers))
    random_ase_uids = sample(ase_data$UID, nrow(ase_outliers))
    risk_df = data.frame(Feature = character(), MedZ_Risk = numeric(), Splice_Risk = numeric(), ASE_Risk = numeric())
    for (vcat in unique(variant_data$variant_cat)) {
        all_with = filter(variant_data,variant_cat==vcat, gene_id %in% medz_data$Gene, indiv_id %in% medz_data$Ind)
        all_with_uids = unique(all_with$UID)
        numYes = length(which(all_with_uids %in% random_medz_uids))
        numNo = length(which(!(all_with_uids %in% random_medz_uids)))
        mrisk = numYes / length(all_with_uids)
        all_with = filter(variant_data,variant_cat==vcat,gene_id %in% splicing_data$CLUSTER_ID, indiv_id %in% splicing_data$variable)
        all_with_uids = unique(all_with$UID)
        numYes = length(which(all_with_uids %in% random_sp_uids))
        numNo = length(which(!(all_with_uids %in% random_sp_uids)))
        srisk = numYes / length(all_with_uids)
        all_with = filter(variant_data,variant_cat==vcat, gene_id %in% ase_data$Gene, indiv_id %in% ase_data$SampleName)
        all_with_uids = unique(all_with$UID)
        numYes = length(which(all_with_uids %in% random_ase_uids))
        numNo = length(which(!(all_with_uids %in% random_ase_uids)))
        arisk = numYes / length(all_with_uids)
        all_risks = data.frame(Feature=vcat,MedZ_Risk=mrisk,Splice_Risk=srisk,ASE_Risk=arisk)
        write.table(all_risks,file=paste0(data_dir, 'absolute_outlier_risk_exp_splice_ase_sv7_allGenes_allVariants_random.txt'),sep='\t',quote=F,row.names=F,append=T)
    }
    return(risk_df)
}

