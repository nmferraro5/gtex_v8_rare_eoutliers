library(data.table)
library(dplyr)
library(ggplot2)
library(epitools)
require(foreach)
require(doMC)

#data_dir = '/Users/nicoleferraro/durga_local/data/goats_data/v8_data/'
data_dir = '/users/nferraro/data/goats_data/v8_data/'
registerDoMC(cores = 5)

eqtl_coloc_results = fread(paste0(data_dir, 'coloc/coloc_sig_results_0.5.txt'))
eqtl_coloc_results$Trait = sapply(eqtl_coloc_results$filename, function(x) strsplit(x, '__PM__')[[1]][1])

eqtl_results = fread(paste0(data_dir, 'coloc/enloc_sig_results_0.5.txt'))
eqtl_results$Trait = sapply(eqtl_results$filename, function(x) strsplit(x, '_w_')[[1]][2])
eqtl_results$Trait = sapply(eqtl_results$Trait, function(x) strsplit(x, '_enloc_output.txt')[[1]][1])
sqtl_results = fread(paste0(data_dir, 'splicing/sqtl_enloc_rcp0.5_collapsed_bygene.txt'))
colnames(sqtl_results) = c('Chr', 'Start', 'End', 'rcp', 'Trait', 'Tissue', 'GeneChr', 'GeneStart', 'GeneEnd', 'Gene')
sqtl_results = sqtl_results %>% select(rcp, Trait, Gene)

te_variants = fread(paste0(data_dir, 'gtexV8.outliers.v8ciseQTLs.globalOutliers.removed.medz.withraresnps.txt'))
ase_variants = fread(paste0(data_dir, 'ASE/gtexV8.outliers.ase.withraresnps.v8.update.txt'))
sp_variants = fread(paste0(data_dir, 'splicing/gtexV8.outliers.splicing.withraresnps.txt'))
all_variants = unique(c(te_variants$VID, ase_variants$VID, sp_variants$VID))

### UKBB phase 2 stats, filtered to all GTEx rare variants
ukbb_sum_stats = fread(paste0(data_dir, 'coloc/gwas/gtex.coloc.traits.beta.gwas.overlap.txt'))
ukbb_sum_stats$Chr = sapply(ukbb_sum_stats$variant, function(x) strsplit(x, ':')[[1]][1])
ukbb_sum_stats$Pos = sapply(ukbb_sum_stats$variant, function(x) strsplit(x, ':')[[1]][2])
ukbb_sum_stats$VID = paste(paste0('chr',ukbb_sum_stats$Chr), ukbb_sum_stats$Pos, sep=':')

trait_table = fread(paste0(data_dir, 'coloc/gwas/ukbb.trait.table.txt'))
colnames(trait_table)[1] = 'Trait'
ukbb_sum_stats = inner_join(ukbb_sum_stats, trait_table, by='Trait')
ukbb_sum_stats = filter(ukbb_sum_stats, !is.na(beta))
ukbb_sum_stats = ukbb_sum_stats %>% group_by(Trait) %>% mutate(VRank = rank(abs(beta))) %>% 
  mutate(VPercentile = ntile(abs(beta),100)) %>% ungroup()

traits_keep = fread(paste0(data_dir, 'coloc/gwas/only.ukbb.traits.coloc.enloc.sqtl.map.txt'))$Trait

ukbb_sum_coloc_stats = ukbb_sum_stats %>% filter(Trait %in% traits_keep) %>%
  mutate(RareTE = ifelse(VID %in% te_variants$VID, 1, 0)) %>%
  mutate(RareASE = ifelse(VID %in% ase_variants$VID, 1, 0)) %>%
  mutate(RareSplice = ifelse(VID %in% sp_variants$VID, 1, 0)) %>%
  mutate(RareAny = ifelse(VID %in% all_variants, 1, 0))

rare_vid_map = fread(paste0(data_dir, 'coloc/gwas/ukbb_rare_SNPs_genes_map.txt'))
ukbb_sum_coloc_stats = merge(ukbb_sum_coloc_stats, rare_vid_map, by='VID', all.x=T)

ukbb.coloc.map = fread(paste0(data_dir, 'coloc/gwas/only.ukbb.traits.coloc.enloc.sqtl.map.txt'))
ukbb_coloc = filter(ukbb_sum_coloc_stats, Trait %in% ukbb.coloc.map$Trait)
ukbb_coloc = inner_join(ukbb_coloc, ukbb.coloc.map, by=c('Trait', 'Description'))

eqtl_results = eqtl_results %>% select(rcp, Trait, gene_id)
colnames(eqtl_results) = c('ColocProb', 'Enloc_name', 'Gene')
eqtl_coloc_results = eqtl_coloc_results %>% select(p4,Trait,gene)
colnames(eqtl_coloc_results) = c('ColocProb', 'Coloc_name', 'Gene')
colnames(sqtl_results) = c('ColocProb', 'Sqtl_name', 'Gene')

ukbb_eqtl_coloc = inner_join(ukbb_coloc, eqtl_coloc_results, by=c('Coloc_name', 'Gene')) %>%
  mutate(QTL = 'Expression')
ukbb_eqtl_enloc = inner_join(ukbb_coloc, eqtl_results, by=c('Enloc_name', 'Gene')) %>%
  mutate(QTL = 'Expression')
ukbb_sqtl = inner_join(ukbb_coloc, sqtl_results, by=c('Sqtl_name', 'Gene')) %>%
  mutate(QTL = 'Splicing')
ukbb_coloc_all = rbind(ukbb_eqtl_coloc, ukbb_eqtl_enloc, ukbb_sqtl) %>%
  group_by(VID,Gene,Trait) %>% sample_n(1) %>% ungroup()

ukbb_sum = ukbb_coloc_all %>% filter(low_confidence_variant == 'FALSE') %>%
  group_by(VID,Trait) %>% sample_n(1) %>% ungroup() %>% 
  select(VID,Gene,Trait,beta,VRank,RareTE,RareASE,RareSplice,RareAny)

## split into outlier, high WS, high CADD categories
exp_ws = fread(paste0(data_dir, 'coloc/gwas/expression.watershed.ukbb.variants.v8.update.txt'))
exp_ws$Gene = sapply(exp_ws$sample_names, function(x) strsplit(x, ':')[[1]][2])
sp_ws = fread(paste0(data_dir, 'coloc/gwas/splicing.watershed.ukbb.variants.v8.update.txt'))
sp_ws$Gene = sapply(sp_ws$sample_names, function(x) strsplit(x, ':')[[1]][2])
ase_ws = fread(paste0(data_dir, 'coloc/gwas/ase.watershed.ukbb.variants.v8.update.txt'))
ase_ws$Gene = sapply(ase_ws$sample_names, function(x) strsplit(x, ':')[[1]][2])

exp_ws = inner_join(exp_ws, ukbb_sum_coloc_stats %>% select(minor_AF,low_confidence_variant,beta,pval,VPercentile,Trait,Description,VID), by='VID')
sp_ws = inner_join(sp_ws, ukbb_sum_coloc_stats %>% select(minor_AF,low_confidence_variant,beta,pval,VPercentile,Trait,Description,VID), by='VID')
ase_ws = inner_join(ase_ws, ukbb_sum_coloc_stats %>% select(minor_AF,low_confidence_variant,beta,pval,VPercentile,Trait,Description,VID), by='VID')

ws_all = rbind(exp_ws %>% mutate(ws_posterior = total_expression_watershed_posterior) %>% 
                 select(-total_expression_watershed_posterior) %>% mutate(Type = 'TE'),
               sp_ws %>% mutate(ws_posterior = splicing_watershed_posterior) %>% 
                 select(-splicing_watershed_posterior) %>% mutate(Type = 'Splicing'),
               ase_ws %>% mutate(ws_posterior = ase_watershed_posterior) %>% 
                 select(-ase_watershed_posterior) %>% mutate(Type = 'ASE'))

cadd_scores = fread(paste0(data_dir, 'coloc/gwas/gtexV8.ukbb.rare.cadd.scores.txt'))
cadd_scores = cadd_scores %>% group_by(VID) %>% top_n(1,RawScore) %>% sample_n(1) %>% ungroup()

high_ws_variants = filter(ws_all, ws_posterior > 0.5) %>% group_by(VID) %>% top_n(1,ws_posterior) %>%
  sample_n(1) %>% ungroup() %>% select(VID,ws_posterior,Gene)
high_cadd_variants = cadd_scores %>% top_n(length(unique(high_ws_variants$VID)), RawScore) %>% select(VID,RawScore)
#high_cadd_variants = cadd_scores %>% top_n(730, RawScore) %>% select(VID)

ukbb_coloc_highconf = filter(ukbb_coloc_all, low_confidence_variant == 'FALSE') %>%
  group_by(Trait) %>% mutate(VPercentile = ntile(abs(beta),100)) %>% ungroup()

ukbb_coloc_high_ws = inner_join(ukbb_coloc_highconf, high_ws_variants, by=c('VID','Gene')) %>%
  group_by(VID,Trait) %>% sample_n(1) %>% ungroup()
bg_coloc_high_ws = filter(ukbb_coloc_highconf, !(VID %in% ukbb_coloc_high_ws$VID),
                          Trait %in% ukbb_coloc_high_ws$Trait) %>% group_by(VID,Trait) %>%
  sample_n(1) %>% ungroup()
ukbb_coloc_high_cadd = filter(ukbb_coloc_highconf,VID %in% high_cadd_variants$VID) %>% group_by(VID,Trait) %>%
  sample_n(1) %>% ungroup()
bg_coloc_high_cadd = filter(ukbb_coloc_highconf, !(VID %in% ukbb_coloc_high_ws$VID),
                            Trait %in% ukbb_coloc_high_cadd$Trait) %>% group_by(VID,Trait) %>%
  sample_n(1) %>% ungroup()


## test max effect size percentile vs Watershed and CADD
ukbb_top_summary = ukbb_sum_stats %>% group_by(VID) %>% top_n(1,VPercentile) %>% top_n(1,-minor_AF) %>% sample_n(1) %>% ungroup()
ws_summary = ws_all %>% group_by(VID,Gene) %>% mutate(MaxWS = max(ws_posterior)) %>%
  sample_n(1) %>% ungroup() %>% select(VID,Gene,MaxWS)
gam_summary = fread(paste0(data_dir, 'coloc/gwas/all.gam.ukbb.variants.txt')) %>%
  group_by(VID,Gene) %>% mutate(MaxGAM = max(gam_posterior,na.rm=T)) %>% sample_n(1) %>% ungroup()

ukbb_coloc_summary = inner_join(ukbb_coloc_highconf, ws_summary, by=c('VID','Gene'))
ukbb_coloc_summary = inner_join(ukbb_coloc_summary, cadd_scores %>% select(RawScore,VID), by='VID')
ukbb_coloc_summary = inner_join(ukbb_coloc_summary, gam_summary, by=c('VID','Gene'))
ukbb_coloc_summary = ukbb_coloc_summary %>% group_by(Trait) %>%
  mutate(VPercentile = ntile(abs(beta), 100)) %>% ungroup()

ukbb_coloc_summary$NormWS = scale(ukbb_coloc_summary$MaxWS)
ukbb_coloc_summary$NormCADD = scale(ukbb_coloc_summary$RawScore)

af_stats = summary(lm(VPercentile ~ minor_AF, data=ukbb_coloc_summary))$coefficients[2,]
ws_stats = summary(lm(VPercentile ~ NormWS, data=ukbb_coloc_summary))$coefficients[2,]
gam_stats = summary(lm(VPercentile ~ MaxGAM, data=ukbb_coloc_summary))$coefficients[2,]
cadd_stats = summary(lm(VPercentile ~ NormCADD, data=ukbb_coloc_summary))$coefficients[2,]
all_stats = as.data.frame(rbind(af_stats, ws_stats, gam_stats, cadd_stats))
all_stats$Feature = c('MAF', 'Watershed', 'GAM', 'CADD')
colnames(all_stats)[2] = 'SE'
all_stats$Feature = factor(all_stats$Feature, levels=c('CADD', 'GAM', 'Watershed', 'MAF'))
lm_plot = ggplot(all_stats %>% filter(Feature != 'MAF', Feature != 'GAM'), aes(x=Feature, y= abs(Estimate))) + geom_point(size=2) + 
  geom_errorbar(aes(ymin = abs(Estimate) - 1.95*SE, ymax = abs(Estimate) + 1.95*SE), width=0) + 
  theme_bw() + geom_hline(yintercept=0, color='grey') + xlab('') + ylab('') +
  ylab('Coloc region variant effect size percentile ~ Feature') +
  annotate("text", x=1, y=1.75, label='*', cex=5) +
  annotate("text", x=2, y=2.5, label='***', cex=5) +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16)) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

mc_plot = ggplot(ukbb_coloc_summary, aes(x=MaxWS, y=RawScore)) + geom_point() + theme_bw() +
  xlab('Max Watershed posterior') + ylab('CADD score') +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16)) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

write.table(ukbb_coloc_summary, file=paste0(data_dir, 'coloc/gwas/ukbb.coloc.watershed.cadd.summary.stats.v8.update.txt'), sep='\t', quote=F, row.names=F)

### Look at proportion to follow up on after filtering
ukbb_coloc_high_ws = inner_join(ukbb_coloc_highconf, high_ws_variants, by=c('VID','Gene')) %>%
  group_by(VID,Gene,Trait) %>% sample_n(1) %>% ungroup()

trait_table = data.frame(Trait = character(), NG = numeric(), NV = numeric(), BestRank = numeric(), BestPercentile = numeric(), BestWS = numeric())
for (gtrait in unique(ukbb_coloc_high_ws$Trait)) {
  tdata = filter(ukbb_coloc_high_ws, Trait == gtrait)
  tname = tdata$Description[1]
  ng = length(unique(tdata$Gene))
  nv = length(unique(tdata$VID))
  brank = max(tdata$VRank)
  bp = max(tdata$VPercentile)
  best_vid = filter(tdata, VRank == max(tdata$VRank))$VID
  best_gene = filter(tdata, VRank == max(tdata$VRank))$Gene
  best_ws = filter(tdata, VRank == max(tdata$VRank))$ws_posterior
  trait_table = rbind(trait_table, data.frame(Trait = tname, NG = ng, NV = nv, BestRank = brank, BestPercentile = bp, BestWS = best_ws))
}

ws_map = ws_all %>% group_by(VID,Gene) %>% top_n(1,ws_posterior) %>% sample_n(1) %>% ungroup()
ukbb_coloc_highconf = inner_join(ukbb_coloc_highconf, ws_map %>% select(VID,Gene,ws_posterior), by=c('VID','Gene'))
ukbb_coloc_highconf = inner_join(ukbb_coloc_highconf, cadd_scores %>% select(VID,RawScore), by=c('VID'))

ws_cadd_plot = ggplot(ukbb_coloc_highconf, aes(x=RawScore,y=ws_posterior)) + geom_point() +
  xlab('CADD Score') + ylab('Watershed posterior') +
  geom_vline(xintercept=2.3,color='firebrick4') + 
  geom_hline(yintercept=0.5,color='firebrick4') +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16)) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ukbb_coloc_highconf = ukbb_coloc_highconf %>% group_by(Trait) %>% mutate(VPercentile = ntile(abs(beta),100)) %>% ungroup()
hit_table = data.frame(PT = numeric(), NumWS = numeric(), NumWSC = numeric(), NumVGT = numeric(), NumWSH = numeric(),
                       CT = numeric(), NumCADDC = numeric(), NumCVGT = numeric(), NumCADDH = numeric(),
                       NumGAMC = numeric(), NumGVGT = numeric(), NumGAMH = numeric())
for (pt in c(0.01, 0.05,0.25,0.5,0.75,0.9)) {
  high_ws_variants = filter(ws_all, ws_posterior > pt) %>% group_by(VID,Gene) %>% top_n(1,ws_posterior) %>%
    sample_n(1) %>% ungroup() %>% select(VID,ws_posterior,Gene)
  high_cadd_variants = cadd_scores %>% top_n(length(unique(high_ws_variants$VID)), RawScore)
  high_gam_variants = gam_summary %>% top_n(length(unique(high_ws_variants$VID)), MaxGAM)
  ct = min(high_cadd_variants$RawScore)
  nws = length(unique(high_ws_variants$VID))
  ukbb_coloc_high_ws = inner_join(ukbb_coloc_highconf, high_ws_variants, by=c('VID','Gene')) %>%
   group_by(VID,Gene,Trait) %>% sample_n(1) %>% ungroup()
  ukbb_coloc_high_gam = inner_join(ukbb_coloc_highconf, high_gam_variants, by=c('VID','Gene')) %>%
    group_by(VID,Gene,Trait) %>% sample_n(1) %>% ungroup()
  # ukbb_coloc_high_ws = inner_join(ukbb_sum_stats, high_ws_variants, by='VID') %>%
  #   group_by(VID,Trait) %>% sample_n(1) %>% ungroup()
  nwsc = length(unique(ukbb_coloc_high_ws$VID))
  ngamc = length(unique(ukbb_coloc_high_gam$VID))
  ukbb_coloc_high_cadd = filter(ukbb_coloc_highconf,VID %in% high_cadd_variants$VID) %>% group_by(VID,Gene,Trait) %>%
    sample_n(1) %>% ungroup()
  # ukbb_coloc_high_cadd = filter(ukbb_sum_stats,VID %in% high_cadd_variants$VID) %>% group_by(VID,Trait) %>%
  #   sample_n(1) %>% ungroup()
  ncaddc = length(unique(ukbb_coloc_high_cadd$VID))
  gam_hp = nrow(ukbb_coloc_high_gam)
  ws_hp = nrow(ukbb_coloc_high_ws)
  cadd_hp = nrow(ukbb_coloc_high_cadd)
  ws_hits = nrow(filter(ukbb_coloc_high_ws, VPercentile >= 75))
  gam_hits = nrow(filter(ukbb_coloc_high_gam, VPercentile >= 75))
  #ws_hits = length(unique((filter(ukbb_coloc_high_ws, VPercentile >= 75))$VID))
  cadd_hits = nrow(filter(ukbb_coloc_high_cadd, VPercentile >= 75))
  #cadd_hits = length(unique((filter(ukbb_coloc_high_cadd, VPercentile >= 75))$VID))
  hit_table = rbind(hit_table, data.frame(PT = pt, NumWS = nws, NumWSC = nwsc, NumVGT = ws_hp, NumWSH = ws_hits,
                         CT = ct, NumCADDC = ncaddc, NumCVGT = cadd_hp, NumCADDH = cadd_hits,
                         NumGAMC = ngamc, NumGVGT = gam_hp, NumGAMH = gam_hits))
}

hit_random_table = data.frame(PT = numeric(), NumWS = numeric(), NumWSC = numeric(), NumVGT = numeric(), NumWSH = numeric(), Model = character())
hit_random_table = foreach(i = 1:1000, .combine = rbind) %dopar% {
  print(i)
  hit_new_random_table = data.frame(PT = numeric(), NumWS = numeric(), NumWSC = numeric(), NumVGT = numeric(), NumWSH = numeric(), Model = character())
  for (pt in c(0.01, 0.05,0.25,0.5,0.75,0.9)) {
    high_ws_variants = filter(ws_all, ws_posterior > pt) %>% group_by(VID,Gene) %>% top_n(1,ws_posterior) %>%
      sample_n(1) %>% ungroup() %>% select(VID,ws_posterior,Gene)
    high_cadd_variants = cadd_scores %>% top_n(length(unique(high_ws_variants$VID)), RawScore)
    high_gam_variants = gam_summary %>% top_n(length(unique(high_ws_variants$VID)), MaxGAM)
    random_variants = sample(unique(ws_all$VID), length(unique(high_ws_variants$VID)))
    random_cadd_variants = sample(unique(ws_all$VID), length(unique(high_cadd_variants$VID)))
    random_gam_variants = sample(unique(gam_summary$VID), length(unique(high_gam_variants$VID)))
    high_ws_variants = filter(ws_all, VID %in% random_variants) %>% group_by(VID,Gene) %>% top_n(1,ws_posterior) %>%
      sample_n(1) %>% ungroup() %>% select(VID,ws_posterior,Gene)
    high_cadd_variants = filter(ws_all, VID %in% random_cadd_variants) %>% group_by(VID,Gene) %>% top_n(1,ws_posterior) %>%
      sample_n(1) %>% ungroup() %>% select(VID,ws_posterior,Gene)
    high_gam_variants = filter(ws_all, VID %in% random_gam_variants) %>% group_by(VID,Gene) %>% top_n(1,ws_posterior) %>%
      sample_n(1) %>% ungroup() %>% select(VID,ws_posterior,Gene)
    nws = length(unique(high_ws_variants$VID))
    ncadd = length(unique(high_cadd_variants$VID))
    ngam = length(unique(high_gam_variants$VID))
    ukbb_coloc_high_ws = inner_join(ukbb_coloc_highconf, high_ws_variants, by=c('VID','Gene')) %>%
      group_by(VID,Gene,Trait) %>% sample_n(1) %>% ungroup()
    ukbb_coloc_high_cadd = inner_join(ukbb_coloc_highconf, high_cadd_variants, by=c('VID','Gene')) %>%
      group_by(VID,Gene,Trait) %>% sample_n(1) %>% ungroup()
    ukbb_coloc_high_gam = inner_join(ukbb_coloc_highconf, high_gam_variants, by=c('VID','Gene')) %>%
      group_by(VID,Gene,Trait) %>% sample_n(1) %>% ungroup()
    # ukbb_coloc_high_ws = inner_join(ukbb_sum_stats, high_ws_variants, by='VID') %>%
    #   group_by(VID,Trait) %>% sample_n(1) %>% ungroup()
    # ukbb_coloc_high_cadd = inner_join(ukbb_sum_stats, high_cadd_variants, by='VID') %>%
    #   group_by(VID,Trait) %>% sample_n(1) %>% ungroup()
    nwsc = length(unique(ukbb_coloc_high_ws$VID))
    ws_hp = nrow(ukbb_coloc_high_ws)
    ws_hits = nrow(filter(ukbb_coloc_high_ws, VPercentile >= 75))
    
    ncsc = length(unique(ukbb_coloc_high_cadd$VID))
    cadd_hp = nrow(ukbb_coloc_high_cadd)
    cadd_hits = nrow(filter(ukbb_coloc_high_cadd, VPercentile >= 75))
    ngsc = length(unique(ukbb_coloc_high_gam$VID))
    gam_hp = nrow(ukbb_coloc_high_gam)
    gam_hits = nrow(filter(ukbb_coloc_high_gam, VPercentile >= 75))
    #ws_hits = length(unique((filter(ukbb_coloc_high_ws, VPercentile >= 75))$VID))
    #cadd_hits = length(unique((filter(ukbb_coloc_high_cadd, VPercentile >= 75))$VID))
    hit_new_random_table = rbind(hit_new_random_table, data.frame(PT = pt, NumWS = nws, NumWSC = nwsc, NumVGT = ws_hp, NumWSH = ws_hits, Model = 'Watershed'))
    hit_new_random_table = rbind(hit_new_random_table, data.frame(PT = pt, NumWS = ncadd, NumWSC = ncsc, NumVGT = cadd_hp, NumWSH = cadd_hits, Model = 'CADD'))
    hit_new_random_table = rbind(hit_new_random_table, data.frame(PT = pt, NumWS = ngam, NumWSC = ngsc, NumVGT = gam_hp, NumWSH = gam_hits, Model = 'GAM'))
  }
  hit_new_random_table
}

hit_ws_table = hit_table[,1:5] %>% mutate(Cat = 'Actual') %>% mutate(Prop = NumWSH/NumVGT)
#hit_ws_table = hit_table[,1:5] %>% mutate(Cat = 'Actual') %>% mutate(Prop = NumWSH/NumWSC)
hit_cadd_table = hit_table[,c(1,7:9)] %>% mutate(Model = 'CADD', Cat = 'Actual') %>% mutate(Prop = NumCADDH/NumCVGT)
hit_gam_table = hit_table[,c(1,10:12)] %>% mutate(Model = 'GAM', Cat = 'Actual') %>% mutate(Prop = NumGAMH/NumGVGT)

#hit_cadd_table = hit_table[,c(1,7:9)] %>% mutate(Model = 'CADD', Cat = 'Actual') %>% mutate(Prop = NumCADDH/NumCADDC)
hit_cadd_random_table = rbind(hit_random_table %>% filter(Model == 'CADD') %>%
                                mutate(Cat = 'Random',Prop = NumWSH/NumVGT) %>% select(PT,Cat,Prop), hit_cadd_table %>% select(PT,Cat,Prop))
#hit_cadd_random_table = rbind(hit_random_table %>% filter(Model == 'CADD') %>%
#                                mutate(Cat = 'Random',Prop = NumWSH/NumWSC) %>% select(PT,Cat,Prop), hit_cadd_table %>% select(PT,Cat,Prop))

hit_gam_random_table = rbind(hit_random_table %>% filter(Model == 'GAM') %>%
                                mutate(Cat = 'Random',Prop = NumWSH/NumVGT) %>% select(PT,Cat,Prop), hit_gam_table %>% select(PT,Cat,Prop))

hit_cadd_random_table$PT = factor(hit_cadd_random_table$PT, levels=c(0.01,0.05,0.25,0.5,0.75,0.9))
hit_ws_random_table = rbind(hit_random_table %>% filter(Model == 'Watershed') %>%
                                 mutate(Cat = 'Random',Prop = NumWSH/NumVGT) %>% select(PT,Cat,Prop), hit_ws_table %>% select(PT,Cat,Prop))
#hit_ws_random_table = rbind(hit_random_table %>% filter(Model == 'Watershed') %>%
#                              mutate(Cat = 'Random',Prop = NumWSH/NumWSC) %>% select(PT,Cat,Prop), hit_ws_table %>% select(PT,Cat,Prop))
hit_ws_random_table$PT = factor(hit_ws_random_table$PT, levels=c(0.01,0.05,0.25,0.5,0.75,0.9))
# ggplot(hit_gam_random_table %>% filter(Cat == 'Random'), aes(x=PT,y=Prop)) + geom_boxplot() +
#   geom_point(data = filter(hit_gam_random_table, Cat == 'Actual'), aes(x=PT, y=Prop),color = 'red',size=2) +
#   theme_bw() + xlab('Posterior threshold') + ylab('Proportion of coloc hits') +
#   theme(axis.text=element_text(size=14),
#         axis.title=element_text(size=16)) +
#   theme(panel.border = element_blank(),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"))

hit_both_random_table = rbind(hit_ws_random_table %>% mutate(Model = 'Watershed'), hit_cadd_random_table %>% mutate(Model = 'CADD'))
#hit_both_random_table$PT = factor(hit_both_random_table$PT, levels=c(0.01,0.05,0.25,0.5,0.75,0.9))
htable = ggplot(hit_both_random_table %>% filter(Cat == 'Random'), aes(x=PT,y=Prop)) + geom_boxplot() +
  geom_point(data = filter(hit_both_random_table, Cat == 'Actual'), aes(x=PT, y=Prop),color = 'red',size=2) +
  theme_bw() + xlab('Posterior threshold') + ylab('Proportion of variants in top 25% of effect sizes') +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16)) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  facet_grid(.~Model)

write.table(hit_both_random_table, file=paste0(data_dir, 'coloc/gwas/gtexV8.ukbb.watershed.proportion.hits.v8.update.txt'), sep='\t', quote=F, row.names=F)


