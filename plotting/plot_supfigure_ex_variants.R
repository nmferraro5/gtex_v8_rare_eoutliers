library(ggplot2)
library(ggthemes)
library(data.table)
library(dplyr)
library(RColorBrewer)
library(cowplot)
library(formattable)

gtex_v8_figure_theme <- function() {
  return(theme(plot.title = element_text(face="plain",size=8), text = element_text(size=8),
               axis.text=element_text(size=7), panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(),panel.background = element_blank(), 
               panel.border = element_blank(),
               axis.line = element_line(colour = "black"), 
               legend.text = element_text(size=7), legend.title = element_text(size=8)))
}


data_dir = '/Users/nicoleferraro/durga_local/data/goats_data/v8_data/'

asthma_vars = fread(paste0(data_dir, 'coloc/gwas/traits/asthma_example_variant.txt'))
bmi_vars = fread(paste0(data_dir, 'coloc/gwas/traits/bmi_example_variant.txt'))
asthma_vars = asthma_vars %>% select(Gene,Ind,Chr,Pos,variant_cat,af_gnomad,splicing_pvalue, exp_pvalue, ase_pvalue, splicing_ws, exp_ws, ase_ws)
bmi_vars = bmi_vars %>% select(Gene,Ind,Chr,Pos,variant_cat,af_gnomad,splicing_pvalue, exp_pvalue, ase_pvalue, splicing_ws, exp_ws, ase_ws)
asthma_vars[,c(6:12)] = round(asthma_vars[,c(6:12)],4)
formattable(asthma_vars %>% select(-af_gnomad))
bmi_vars[,c(6:12)] = round(bmi_vars[,c(6:12)],4)
formattable(bmi_vars %>% select(-af_gnomad) %>% filter(Pos == 846407))

## save figure
fig_1 = plot_grid(sfig_A, sfig_B, nrow = 2, labels=c('A', 'B'), align='v')

ggsave(fig_1, file=paste0(data_dir, 'paper_figures/sfig_enrich_v2.pdf'), width=7.2, height=7.2,units="in")


testing = fread(paste0(data_dir, 'gtexV8.all.data.types.snps.indels.continuous.10kb.ase.v8.vg.gnomad.txt'))
testing$Maf = sapply(testing$Maf, function(x)
  ifelse(x == 'MAFnovel', 'novel',
         ifelse(x == 'MAFac1', 'single',
                ifelse(x == 'MAFac2', 'double',
                       ifelse(x == 'MAFgnomad0-1only', 'rare', 'low frequency')))))
testing = filter(testing, Maf != 'single', Maf != 'double')
testing$Maf = factor(testing$Maf, levels=c('novel', 'rare', 'low frequency'))
ggplot(testing, aes(x=Maf, y=Beta)) + geom_point(size=3) + theme_bw() +
  geom_errorbar(aes(ymin = Beta-1.95*SE, ymax = Beta+1.95*SE), width=0) +
  geom_hline(yintercept=0, color='grey') + xlab('') +
  ggtitle('ASE - SNPs and indels') + theme(axis.text=element_text(size=14),
                         axis.title=element_text(size=16),
                         title=element_text(size=16))

risks = c(1.002014, 1.118507, 1.217059, 1.282187)
lowers = c(0.9986417, 1.102581, 1.185706, 1.234683)
uppers = c(1.005398, 1.134664, 1.249241, 1.331518)
risk_df = data.frame(Risks = risks, Lower = lowers, Upper = uppers, Z = c(1,2,3,4))
ggplot(risk_df, aes(x=Z,y=Risks)) +
  geom_point(size=3) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width=0) +
  theme_bw() + ylab('Relative risk') + xlab('|Median Z| threshold') +
  geom_hline(yintercept=1, linetype="dashed", colour="darkgrey") +
  theme(axis.text=element_text(size=14),
        axis.title=element_text(size=16))






