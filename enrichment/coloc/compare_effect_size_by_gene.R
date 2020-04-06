library(data.table)
library(dplyr)
library(ggplot2)
require(doMC)
require(foreach)

registerDoMC(cores = 5)

data_dir = '/users/nferraro/data/goats_data/v8_data/coloc/gwas/'

gtex_gwas = fread(paste0(data_dir, 'gtex.traits.beta.gwas.overlap.txt'))
gtex_gwas$Chr = sapply(gtex_gwas$variant, function(x) strsplit(x, ':')[[1]][1])
gtex_gwas$Pos = sapply(gtex_gwas$variant, function(x) strsplit(x, ':')[[1]][2])
gtex_gwas$VID = paste(paste0('chr',gtex_gwas$Chr), gtex_gwas$Pos, sep=':')
trait_table = fread(paste0(data_dir, 'ukbb.trait.table.txt'))
colnames(trait_table)[1] = 'Trait'
gtex_gwas = inner_join(gtex_gwas, trait_table, by='Trait')

exp_ws = fread(paste0(data_dir, 'expression.watershed.ukbb.variants.txt'))
sp_ws = fread(paste0(data_dir, 'splicing.watershed.ukbb.variants.txt'))
ase_ws = fread(paste0(data_dir, 'ase.watershed.ukbb.variants.txt'))

exp_ws = inner_join(exp_ws, gtex_gwas %>% select(minor_AF,low_confidence_variant,beta,pval,Trait,Description,VID), by='VID')
sp_ws = inner_join(sp_ws, gtex_gwas %>% select(minor_AF,low_confidence_variant,beta,pval,Trait,Description,VID), by='VID')
ase_ws = inner_join(ase_ws, gtex_gwas %>% select(minor_AF,low_confidence_variant,beta,pval,Trait,Description,VID), by='VID')

rm(gtex_gwas)

ws_all = rbind(exp_ws %>% mutate(posterior = total_expression_watershed_posterior) %>% select(-total_expression_watershed_posterior) %>% mutate(Type = 'TE'),
               sp_ws %>% mutate(posterior = splicing_watershed_posterior) %>% select(-splicing_watershed_posterior) %>% mutate(Type = 'Splicing'),
               ase_ws %>% mutate(posterior = ase_watershed_posterior) %>% select(-ase_watershed_posterior) %>% mutate(Type = 'ASE'))

rm(exp_ws, sp_ws, ase_ws)

ws_all$Gene = sapply(ws_all$sample_names, function(x) strsplit(x, ':')[[1]][2])
ws_ranks = ws_all %>% filter(!is.na(beta)) %>% group_by(Trait,VID) %>% 
  sample_n(1) %>% ungroup() %>%
  group_by(Trait) %>% mutate(TRank = rank(abs(beta))) %>% ungroup() %>%
  select(Trait,VID,TRank)

ws_all = inner_join(ws_all, ws_ranks %>% select(Trait,VID,TRank), by=c('Trait','VID'))
trait_table = fread(paste0(data_dir, 'ukbb.gtex.trait.table.txt'))
colnames(trait_table)[1] = 'Trait'
all_pts = c(0.5, 0.6, 0.7, 0.8, 0.9)

get_pt_data <- function(pt) {
  ws_top = filter(ws_all, posterior > pt) %>% select(-sample_names)
  #ws_other = filter(ws_all, Gene %in% ws_top$Gene, !(VID %in% ws_top$VID)) %>% select(-sample_names)
  ws_other = filter(ws_all, !(VID %in% ws_top$VID)) %>% select(-sample_names)
  ws_both = inner_join(ws_top, ws_other, by=c('Gene','Trait', 'Description')) %>%
    mutate(DAF = abs(minor_AF.x - minor_AF.y)) %>% filter(DAF < 0.01) %>%
    select(Gene,VID.x,VID.y,TRank.x,TRank.y,Trait,DAF)

  ws_melted = melt(ws_both, id.vars = c('Gene','VID.x', 'VID.y', 'Trait', 'DAF'))
  ws_melted = inner_join(ws_melted, trait_table %>% select(-Description), by='Trait')

  ws_melted = ws_melted %>% group_by(Gene,VID.x,variable,Name) %>% sample_n(1) %>% ungroup()
  ws_melted$PT = pt
  return(ws_melted)
}

ws_all_melted = do.call(rbind, lapply(all_pts, get_pt_data))
ws_all_melted$PT = factor(ws_all_melted$PT, levels=c(0.5,0.6,0.7,0.8,0.9))
ws_all_melted$variable = sapply(ws_all_melted$variable, function(x)
  ifelse(x == 'TRank.x', 'High Watershed', 'Other'))
ggplot(ws_all_melted %>% filter(PT != 0.5 & PT != 0.6), aes(x=PT, y=value)) + geom_boxplot(aes(fill=variable)) +
  theme_bw() + xlab('Watershed posterior threshold') + ylab('Variant rank') +
  annotate("text", x=1, y=23000, label='***', cex=7) +
  annotate("text", x=2, y=23000, label='***', cex=7) +
  annotate("text", x=3, y=24000, label='****', cex=7) +
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=20),
        legend.text=element_text(size=24),
        legend.title=element_blank()) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))


wilcox.test(filter(ws_all_melted, variable == 'High Watershed', PT == 0.7)$value,
            filter(ws_all_melted, variable == 'Other', PT == 0.7)$value, alternative='g')


## for high watershed variants, get best effect size and see if better than expected by chance
vid_posterior_map = ws_all %>% select(sample_names,posterior,Type) %>% group_by(sample_names,posterior,Type) %>%
  sample_n(1) %>% ungroup()
ws_permute = ws_all %>% select(sample_names, VID, Trait, TRank, beta, posterior, Type)

### permute posteriors within each group
get_random_values <- function(pt) {
  random_values = foreach(i = 1:100, .combine = rbind) %dopar% {
    print(i)
    random_map = vid_posterior_map %>% group_by(Type) %>%
      mutate(NewPosterior = sample(posterior)) %>% ungroup() %>% select(-posterior)
    ws_random = inner_join(ws_permute, random_map, by=c('sample_names', 'Type')) %>% select(-posterior,-Type)
    rm(random_map)
    random_top = filter(ws_random, NewPosterior > pt) %>% group_by(VID) %>% top_n(1,TRank) %>% 
      ungroup() %>% group_by(VID,Trait) %>% top_n(1,NewPosterior) %>% sample_n(1) %>% ungroup()
    colnames(random_top)[6] = 'posterior'
    random_top$PT = pt
    random_top$Iteration = i
    random_top
  }
  return(random_values)
}

ws_top7 = filter(ws_all, posterior > 0.7) %>% group_by(VID) %>% top_n(1,TRank) %>% 
  ungroup() %>% group_by(VID,Trait) %>% top_n(1,posterior) %>% sample_n(1) %>% ungroup()
ws_top8 = filter(ws_all, posterior > 0.8) %>% group_by(VID) %>% top_n(1,TRank) %>% 
  ungroup() %>% group_by(VID,Trait) %>% top_n(1,posterior) %>% sample_n(1) %>% ungroup()
ws_top9 = filter(ws_all, posterior > 0.9) %>% group_by(VID) %>% top_n(1,TRank) %>% 
  ungroup() %>% group_by(VID,Trait) %>% top_n(1,posterior) %>% sample_n(1) %>% ungroup()

actual_values = rbind(ws_top7 %>% mutate(PT = 0.7, Iteration = 0),
                      ws_top8 %>% mutate(PT = 0.8, Iteration = 0),
                      ws_top9 %>% mutate(PT = 0.9, Iteration = 0))

actual_values = actual_values %>% select(sample_names, VID, Trait, TRank, beta, posterior, PT, Iteration)

random_values = rbind(get_random_values(0.7), get_random_values(0.8), get_random_values(0.9))

write.table(rbind(actual_values, random_values), file=paste0(data_dir, 'watershed.summary.actual.random.trait.best.rank.txt'), sep='\t', quote=F, row.names=F)

best_ranks = fread(paste0(data_dir, 'watershed.summary.actual.random.trait.best.rank.txt'))
best_ranks$Cat = sapply(best_ranks$Iteration, function(x) ifelse(x == 0, 'Actual', 'Random'))
best_ranks$PT = factor(best_ranks$PT, levels=c(0.7, 0.8, 0.9))
ggplot(best_ranks, aes(x=PT,y=TRank)) + geom_boxplot(aes(fill=Cat)) + theme_bw()

wilcox.test(filter(best_ranks, PT == 0.9, Cat == 'Actual')$TRank,
            filter(best_ranks, PT == 0.9, Cat == 'Random')$TRank,
            alternative='g')





