library(data.table)
library(dplyr)
library(ggplot2)
# require(doMC)
# require(foreach)
# 
# registerDoMC(cores = 5)

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

exp_ws = inner_join(exp_ws, gtex_gwas %>% select(minor_AF,low_confidence_variant,beta,pval,Trait,Description,VID), by='VID')
sp_ws = inner_join(sp_ws, gtex_gwas %>% select(minor_AF,low_confidence_variant,beta,pval,Trait,Description,VID), by='VID')
ase_ws = inner_join(ase_ws, gtex_gwas %>% select(minor_AF,low_confidence_variant,beta,pval,Trait,Description,VID), by='VID')

ws_all = rbind(exp_ws %>% mutate(outlier_pvalue = total_expression_outlier_pvalue, gam_posterior = total_expression_gam_posterior, ws_posterior = total_expression_watershed_posterior) %>% 
                 select(-total_expression_outlier_pvalue, -total_expression_gam_posterior, -total_expression_watershed_posterior) %>% mutate(Type = 'TE'),
               sp_ws %>% mutate(outlier_pvalue = splicing_outlier_pvalue, gam_posterior = splicing_gam_posterior, ws_posterior = splicing_watershed_posterior) %>% 
                 select(-splicing_outlier_pvalue, -splicing_gam_posterior, -splicing_watershed_posterior) %>% mutate(Type = 'Splicing'),
               ase_ws %>% mutate(outlier_pvalue = ase_outlier_pvalue, gam_posterior = ase_gam_posterior, ws_posterior = ase_watershed_posterior) %>% 
                 select(-ase_outlier_pvalue, -ase_gam_posterior, -ase_watershed_posterior) %>% mutate(Type = 'ASE'))

rm(exp_ws, sp_ws, ase_ws)

all_vids = unique(ws_all$VID)
get_random_medians <- function(NV,gtrait) {
  kvids = sample(all_vids, NV)
  lp_trait_data = filter(ws_all, Trait == gtrait, VID %in% kvids)
  return(median(unique(abs(lp_trait_data$beta)),na.rm=T))
}

pt = 0.9
trait_greater = data.frame(Trait = character(), PropMedian = numeric(), Type = character())
for (gtrait in unique(ws_all$Trait)) {
  print(gtrait)
  print('All')
  hp_trait_data = filter(ws_all, !is.na(outlier_pvalue), Trait == gtrait, ws_posterior > pt) %>%
    group_by(VID) %>% sample_n(1) %>% ungroup()
  random_medians = sapply(1:100, function(x) get_random_medians(nrow(hp_trait_data), gtrait))
  prop = length(which(random_medians > median(abs(hp_trait_data$beta))))
  
  print('Data type specific')
  te_trait_data = filter(ws_all, !is.na(outlier_pvalue), Type == 'TE', Trait == gtrait, ws_posterior > pt) %>%
    group_by(VID) %>% sample_n(1) %>% ungroup()
  random_medians = sapply(1:100, function(x) get_random_medians(nrow(te_trait_data), gtrait))
  te_prop = length(which(random_medians > median(abs(te_trait_data$beta))))
  
  ase_trait_data = filter(ws_all, !is.na(outlier_pvalue), Type == 'ASE', Trait == gtrait, ws_posterior > pt) %>%
    group_by(VID) %>% sample_n(1) %>% ungroup()
  random_medians = sapply(1:100, function(x) get_random_medians(nrow(ase_trait_data), gtrait))
  ase_prop = length(which(random_medians > median(abs(ase_trait_data$beta))))
  
  sp_trait_data = filter(ws_all, !is.na(outlier_pvalue), Type == 'Splicing', Trait == gtrait, ws_posterior > pt) %>%
    group_by(VID) %>% sample_n(1) %>% ungroup()
  random_medians = sapply(1:100, function(x) get_random_medians(nrow(sp_trait_data), gtrait))
  sp_prop = length(which(random_medians > median(abs(sp_trait_data$beta))))
  trait_greater = rbind(trait_greater, data.frame(Trait = gtrait, PropMedian = prop/1000, Type = 'All'),
                        data.frame(Trait = gtrait, PropMedian = te_prop/1000, Type = 'TE'),
                        data.frame(Trait = gtrait, PropMedian = ase_prop/1000, Type = 'ASE'),
                        data.frame(Trait = gtrait, PropMedian = sp_prop/1000, Type = 'Splicing'))
}

write.table(trait_greater, file=paste0(data_dir, 'gtexV8.watershed.ukbb.random.background.variant.effect.size.txt'), sep='\t', quote=F, row.names=F)

trait_greater$Type = factor(trait_greater$Type, levels=c('All', 'TE', 'ASE', 'Splicing'))
ggplot(trait_greater, aes(x=Type,y=PropMedian)) + geom_boxplot() + theme_bw() +
  xlab('') + geom_hline(yintercept=0.05,color='darkgrey', linetype='dashed') +
  ylab('Proportion random samples > actual median per trait') +
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=16)) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

# get_sum_stats <- function(pt) {
#   high_means = c()
#   low_means = c()
#   high_medians = c()
#   low_medians = c()
#   for (gtrait in unique(ws_all$Trait)) {
#     hp_trait_data = filter(ws_all, !is.na(outlier_pvalue), Trait == gtrait, ws_posterior > pt) %>%
#       group_by(VID) %>% sample_n(1) %>% ungroup()
#     lp_trait_data = filter(ws_all, Trait == gtrait, ws_posterior < pt) %>%
#       filter(!(VID %in% hp_trait_data$VID))  %>%
#       group_by(VID) %>% sample_n(1) %>% ungroup()
#     high_means = c(high_means, mean(abs(hp_trait_data$beta),na.rm=T))
#     low_means = c(low_means, mean(abs(lp_trait_data$beta),na.rm=T))
#     high_medians = c(high_medians, median(abs(hp_trait_data$beta),na.rm=T))
#     low_medians = c(low_medians, median(abs(lp_trait_data$beta),na.rm=T))
#   }
# 
#   all_vals = as.data.frame(cbind(high_means, low_means, high_medians, low_medians))
#   all_vals$high_means = as.numeric(as.character(all_vals$high_means))
#   all_vals$low_means = as.numeric(as.character(all_vals$low_means))
#   all_vals$high_medians = as.numeric(as.character(all_vals$high_medians))
#   all_vals$low_medians = as.numeric(as.character(all_vals$low_medians))
#   return(all_vals %>% mutate(PT = pt, Iteration = 0))
# }
# 
# all_pts = c(0.5, 0.7, 0.8, 0.9)
# pt_values = do.call(rbind, lapply(all_pts, get_sum_stats))
# print('Done actual')
# pt_values$Trait = rep(unique(ws_all$Trait), 4)
# 
# vid_posterior_map = ws_all %>% select(sample_names,ws_posterior,Type) %>% group_by(sample_names,ws_posterior,Type) %>%
#   sample_n(1) %>% ungroup()
# 
# ws_permute = ws_all %>% select(sample_names, VID, Trait, beta, ws_posterior, Type)
# rm(ws_all)
# 
# ### permute posteriors within each group
# get_random_values <- function(pt) {
#   random_values = foreach(i = 1:100, .combine = rbind) %dopar% {
#     print(i)
#     random_map = vid_posterior_map %>% group_by(Type) %>%
#       mutate(NewPosterior = sample(ws_posterior)) %>% ungroup() %>% select(-ws_posterior)
#     ws_random = inner_join(ws_permute, random_map, by=c('sample_names', 'Type')) %>% select(-ws_posterior,-Type)
#     rm(random_map)
#     high_means = c()
#     low_means = c()
#     high_medians = c()
#     low_medians = c()
#     for (gtrait in unique(ws_random$Trait)) {
#       hp_trait_data = filter(ws_random, Trait == gtrait, NewPosterior > pt) %>% distinct() %>%
#         select(beta)
#       lp_trait_data = filter(ws_random, NewPosterior < pt) %>%
#         filter(!(VID %in% hp_trait_data$VID)) %>% distinct() %>% select(beta)
#       high_means = c(high_means, mean(abs(hp_trait_data$beta),na.rm=T))
#       low_means = c(low_means, mean(abs(lp_trait_data$beta),na.rm=T))
#       high_medians = c(high_medians, median(abs(hp_trait_data$beta),na.rm=T))
#       low_medians = c(low_medians, median(abs(lp_trait_data$beta),na.rm=T))
#     }
#     all_vals = as.data.frame(cbind(high_means, low_means, high_medians, low_medians))
#     all_vals$high_means = as.numeric(as.character(all_vals$high_means))
#     all_vals$low_means = as.numeric(as.character(all_vals$low_means))
#     all_vals$high_medians = as.numeric(as.character(all_vals$high_medians))
#     all_vals$low_medians = as.numeric(as.character(all_vals$low_medians))
#     all_vals = all_vals %>% mutate(PT = pt, Iteration = i)
#     all_vals$Trait = unique(ws_random$Trait)
#     all_vals
#   }
#   return(random_values)
# }
# 
# random_values7 = get_random_values(0.7)
# write.table(rbind(pt_values %>% mutate(Cat = 'Actual'), random_values7 %>% mutate(Cat = 'Random')), file=paste0(data_dir, 'watershed.summary.actual.random.trait.betas.updated.txt'), sep='\t',quote=F,row.names=F)
# 
# rm(pt_values, random_values7)
# 
# random_values8 = get_random_values(0.8)
# write.table(random_values8 %>% mutate(Cat = 'Random'), file=paste0(data_dir, 'watershed.summary.actual.random.trait.betas.updated.txt'), sep='\t',quote=F,row.names=F,col.names=F, append=T)
# 
# rm(random_values8)
# random_values9 = get_random_values(0.9)
# write.table(random_values9 %>% mutate(Cat = 'Random'), file=paste0(data_dir, 'watershed.summary.actual.random.trait.betas.updated.txt'), sep='\t',quote=F,row.names=F,col.names=F, append=T)
# 
# random_values5 = get_random_values(0.5)
# write.table(random_values5 %>% mutate(Cat = 'Random'), file=paste0(data_dir, 'watershed.summary.actual.random.trait.betas.updated.txt'), sep='\t',quote=F,row.names=F,col.names=F,append=T)
# 
# 
# all_values = fread(paste0(data_dir, 'watershed.summary.actual.random.trait.betas.updated.txt'))
# sum_values = all_values %>% group_by(Iteration,PT) %>%
#   filter(PT != 0.5) %>%
#   mutate(PropMeanHigh = length(which(high_means > low_means))/n()) %>%
#   mutate(PropMedHigh = length(which(high_medians > low_medians))/n()) %>%
#   sample_n(1) %>% ungroup() %>% filter(PT != 'PT')
# 
# sum_values$Cat = sapply(sum_values$Iteration, function(x) ifelse(x == 0, 'Actual', 'Permuted'))
# 
# sum_values$PT = factor(sum_values$PT, levels=c(0.7,0.8,0.9))
# 
# p = ggplot(sum_values, aes(x=PT, y=PropMedHigh)) +
#   geom_jitter(aes(color=Cat),size=4) +
#   scale_color_manual(values=c('darkblue','grey')) +
#   ylim(c(0,1)) +
#   theme_bw() + xlab('Watershed posterior threshold') +
#   ylab('Proportion traits with higher median effect size') +
#   theme(axis.title = element_text(size=22),
#         axis.text = element_text(size=22),
#         legend.text=element_text(size=22),
#         legend.title=element_blank(),
#         legend.position = c(0.1, 0.2)) +
#   theme(panel.border = element_blank(),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"))
# 
# med_vals = melt(all_values %>% select(-high_means, -low_means, -low_medians, -Iteration), id.vars=c('PT','Cat','Trait')) %>%
#   filter(PT == 0.9)
# prop_greater = sapply(unique(med_vals$Trait), function(x)
#   length(which(filter(med_vals, Trait == x, Cat == 'Actual')$value > filter(med_vals, Trait == x, Cat != 'Actual')$value)))
# trait_prop_greater = as.data.frame(cbind(unique(med_vals$Trait), prop_greater))
# trait_table = fread(paste0(data_dir, 'ukbb.gtex.trait.table.txt'))
# colnames(trait_table)[1] = 'Trait'
# colnames(trait_prop_greater) = c('Trait', 'Prop')
# trait_prop_greater = inner_join(trait_prop_greater, trait_table, by='Trait')
# trait_prop_greater$Prop = as.numeric(as.character(trait_prop_greater$Prop)) / 100
# trait_prop_greater = trait_prop_greater %>% arrange(by=Prop)
# trait_prop_greater$Name = factor(trait_prop_greater$Name, levels=unique(trait_prop_greater$Name))
# ggplot(trait_prop_greater, aes(x=Name,y=Prop)) + geom_point() + theme_bw() +
#   ylab('Proportion median effect size greater than random') + xlab('') +
#   geom_hline(yintercept = 0.95, color='grey') +
#   theme(axis.text.x=element_text(size=12,hjust=1,angle=45),
#         axis.text.y=element_text(size=18),
#         axis.title=element_text(size=20)) +
#   theme(panel.border = element_blank(),
#         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#         panel.background = element_blank(), axis.line = element_line(colour = "black"))



