### analysis for supplemental figure for CADD

library(ggplot2)
library(ggthemes)
library(data.table)
library(dplyr)
library(RColorBrewer)
library(cowplot)

gtex_v8_figure_theme <- function() {
  return(theme(plot.title = element_text(face="plain",size=8), text = element_text(size=8),
               axis.text=element_text(size=7), panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(),panel.background = element_blank(), 
               axis.line = element_line(colour = "black"), 
               panel.border=element_blank(),
               legend.text = element_text(size=7), legend.title = element_text(size=8)))
}

data_dir = '/Users/nicoleferraro/durga_local/data/goats_data/v8_data/'

ukbb_coloc_summary = fread(paste0(data_dir, 'coloc/gwas/ukbb.coloc.watershed.cadd.summary.stats.v8.update.txt'))
sfig_A = ggplot(ukbb_coloc_summary, aes(MaxWS)) + geom_histogram() + 
  xlab('Watershed posterior') + ylab('Count') + gtex_v8_figure_theme()
sfig_B = ggplot(ukbb_coloc_summary, aes(RawScore)) + geom_histogram() + 
  xlab('CADD score') + ylab('Count') + gtex_v8_figure_theme()

ukbb_coloc_highconf = fread(paste0(data_dir, 'coloc/gwas/ukbb.coloc.highconf.ws.cadd.v8.update.txt'))
sfig_C = ggplot(ukbb_coloc_highconf, aes(x=RawScore,y=ws_posterior)) + geom_point() +
  xlab('CADD Score') + ylab('Watershed posterior') +
  geom_vline(xintercept=2.3,color='blue') + 
  geom_hline(yintercept=0.5,color='blue') +
  gtex_v8_figure_theme()

ws.cadd.types = fread(paste0(data_dir, 'coloc/gwas/high.ws.cadd.variant.list.vep.txt'))
ws.types = filter(ws.cadd.types, V11 == 'Watershed') %>% distinct()
cadd.types = filter(ws.cadd.types, V11 == 'CADD') %>% distinct()

ws.type.summary = ws.types %>% mutate(N=n()) %>% group_by(V5) %>% mutate(NT = n()) %>% 
  mutate(TypeProp = NT/N) %>% sample_n(1) %>% ungroup()
cadd.type.summary = cadd.types %>% mutate(N=n()) %>% group_by(V5) %>% mutate(NT = n()) %>% 
  mutate(TypeProp = NT/N) %>% sample_n(1) %>% ungroup()
combo.summary = inner_join(ws.type.summary, cadd.type.summary, by='V5') %>%
  mutate(FC = TypeProp.x / TypeProp.y)
combo.summary$Name = c("3' UTR", "downstream", "intergenic", "intronic", 'missense', 'non-coding transcript exon', 'regulatory', 'splice', 'synonymous', 'upstream')
sfig_D = ggplot(combo.summary, aes(x=Name, y=FC)) + geom_bar(stat='identity') + 
  scale_y_log10() + theme_bw() + xlab('') +
  ylab('Proportion watershed / Proportion CADD') +
  gtex_v8_figure_theme() +
  theme(axis.text.x=element_text(hjust=1,angle=45))

hit_cadd_random_table = fread(paste0(data_dir, 'coloc/gwas/gtexV8.ukbb.watershed.proportion.hits.v8.update.txt'))
hit_cadd_random_table$PT = factor(hit_cadd_random_table$PT, levels=c(0.01,0.05,0.25,0.5,0.75,0.9))
sfig_E = ggplot(hit_cadd_random_table %>% filter(Cat == 'Random', Model == 'CADD'), aes(x=PT,y=Prop)) + geom_boxplot() +
  geom_point(data = filter(hit_cadd_random_table, Cat == 'Actual', Model == 'CADD'), aes(x=PT, y=Prop),color = 'red',size=2) +
  theme_bw() + xlab('Posterior threshold') + ylab('Proportion of variants in top 25%') +
  gtex_v8_figure_theme()

first_row <- plot_grid(sfig_A, sfig_B, sfig_C, ncol=3,labels=c('A', 'B', 'C'), align='hv', axis='tlbr')
second_row <- plot_grid(sfig_D, sfig_E, ncol=2, labels=c('D', 'E'), align='hv', axis='tlbr')

sfig <- plot_grid(first_row, second_row, nrow=2)

ggsave(sfig, file=paste0(data_dir, 'paper_figures/revisions/supp/sfig_cadd.png'), width=7.2, height=8, units="in")
