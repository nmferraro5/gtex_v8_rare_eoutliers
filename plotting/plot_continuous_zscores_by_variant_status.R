library(data.table)
library(dplyr)
library(ggpubr)

data_dir = '/users/nferraro/data/goats_data/v8_data/'

data = fread(paste0(data_dir, 'gtexV8.all.outliers.with.rare.variant.status.txt'))

gtex_v8_figure_theme <- function() {
  return(theme(plot.title = element_text(face="plain",size=10), text = element_text(size=10),
               axis.text=element_text(size=9), panel.grid.major = element_blank(), 
               panel.grid.minor = element_blank(),panel.background = element_blank(), 
               axis.line = element_line(colour = "black"), 
               panel.border=element_blank(),
               legend.text = element_text(size=9), legend.title = element_text(size=10)))
}

get_variant_controls <- function(vtype, dtype) {
  print(vtype)
  print(dtype)
  wvariant = filter(data, VarCat == vtype, Cat == dtype) %>% mutate(Status=1)
  wovariant = filter(data, is.na(VarCat), Cat == dtype, Gene %in% wvariant$Gene) %>%
    mutate(Status=0)
  wovariant$VarCat = vtype
  var_data = rbind(wvariant, wovariant)
  return(var_data %>% select(OutlierValue,Cat,VarCat,Status))
}

# vtypes = c('SNPs', 'indels', 'SVs')
# dtypes = c('eOutliers', 'sOutliers', 'aseOutliers')
# types = rbind(rep(vtypes,3), rep(dtypes,each=3))
# all_var_data = do.call(rbind, lapply(1:ncol(types), function(x) get_variant_controls(types[1,x], types[2,x])))
# write.table(all_var_data, file=paste0(data_dir, 'gtexV8.all.outliers.with.rare.variant.status.filtered.txt'), sep='\t', quote=F, row.names=F)

data_dir = '/Users/nicoleferraro/durga_local/data/goats_data/v8_data/'
all_var_data = fread(paste0(data_dir, 'gtexV8.outlier.controls.all.data.types.with.10kb.novel.variant.status.filtered.txt'))
all_var_data = all_var_data %>% mutate(OutlierValue = ifelse(Category == 'eOutliers', -log10(2*pnorm(-abs(OutlierValue))), -log10(OutlierValue)))
all_var_data = filter(all_var_data, !is.na(OutlierValue))
sv_data = filter(all_var_data, InSV == 1) %>% select(-HasNovelSNP, -HasNovelIndel)
all_var_data$HasNovelIndel = factor(all_var_data$HasNovelIndel, levels=c(0,1))
all_var_data$HasNovelSNP = factor(all_var_data$HasNovelSNP, levels=c(0,1))
sv_data$HasNovelSV = factor(sv_data$HasNovelSV, levels=c(0,1))
all_var_data$OutlierValue = sapply(all_var_data$OutlierValue, function(x) ifelse(is.infinite(x), 6.3, x))
sv_data$OutlierValue = sapply(sv_data$OutlierValue, function(x) ifelse(is.infinite(x), 6.3, x))

all_var_data$HasSI = sapply(1:nrow(all_var_data), function(x) ifelse(all_var_data$HasNovelIndel[x] == 1 | all_var_data$HasNovelSNP[x] == 1, 1, 0))

exp_snps = summary(glm(HasSI ~ OutlierValue, data = filter(all_var_data, Category == 'eOutliers'), family = binomial))$coefficients
sp_snps = summary(glm(HasSI ~ OutlierValue, data = filter(all_var_data, Category == 'sOutliers'), family = binomial))$coefficients
ase_snps = summary(glm(HasSI ~ OutlierValue, data = filter(all_var_data, Category == 'aseOutliers'), family = binomial))$coefficients

exp_svs = summary(glm(HasNovelSV ~ OutlierValue, data = filter(sv_data, Category == 'eOutliers'), family = binomial))$coefficients
sp_svs = summary(glm(HasNovelSV ~ OutlierValue, data = filter(sv_data, Category == 'sOutliers'), family = binomial))$coefficients
ase_svs = summary(glm(HasNovelSV ~ OutlierValue, data = filter(sv_data, Category == 'aseOutliers'), family = binomial))$coefficients

var_coefs = rbind(as.data.frame(exp_snps) %>% mutate(Category = 'eOutliers', Variant = 'SNVs+indels'),
                  as.data.frame(sp_snps) %>% mutate(Category = 'sOutliers', Variant = 'SNVs+indels'),
                  as.data.frame(ase_snps) %>% mutate(Category = 'aseOutliers', Variant = 'SNVs+indels'),
                  as.data.frame(exp_svs) %>% mutate(Category = 'eOutliers', Variant = 'SVs'),
                  as.data.frame(sp_svs) %>% mutate(Category = 'sOutliers', Variant = 'SVs'),
                  as.data.frame(ase_svs) %>% mutate(Category = 'aseOutliers', Variant = 'SVs'))

var_coefs = var_coefs[c(2,4,6,8,10,12),]
colnames(var_coefs)[1:2] = c('Beta', 'SE')
var_coefs$Category = factor(var_coefs$Category, levels=c('eOutliers', 'aseOutliers', 'sOutliers'))
sf1 = ggplot(var_coefs, aes(x=Category,y=Beta)) + geom_point(size=2) +
  geom_errorbar(aes(ymin = Beta - 1.95*SE, ymax = Beta + 1.95*SE), width=0,position=position_dodge(width=0.5)) +
  geom_hline(yintercept=0,color='grey',linetype='dashed') + xlab('') +
  theme_bw() + gtex_v8_figure_theme() +
  facet_wrap(.~Variant,scales='free') + theme(strip.background = element_blank())

load(paste0(data_dir, 'gtexV8.continuous.risks.all.types.RData'))
all_coefs$Variant = factor(all_coefs$Variant, levels=c('no_variant','other_noncoding','TE', 'coding','TSS', 'conserved_noncoding','INV','BND','DEL','CNV','DUP', 'splice_region_variant', 'splice_acceptor_variant','frameshift','splice_donor_variant', 'stop'))
plot_cols['no_variant'] = 'darkgrey'
all_coefs$Category = factor(all_coefs$Category, levels=c('eOutliers', 'aseOutliers', 'sOutliers'))

sf2 = ggplot(all_coefs, aes(x=Variant,y=Beta,Group=Category)) +
  geom_point(size=2, aes(color=Variant,shape=Category),position=position_dodge(width=0.5)) +
  scale_shape_manual(values=c('square', 'circle', 'triangle')) +
  geom_errorbar(aes(ymin = Beta - 1.95*SE, ymax = Beta + 1.95*SE), width=0,position=position_dodge(width=0.5)) +
  theme_bw() + ylab('Beta') + xlab('') + ylim(c(-0.5,1.5)) +
  ggtitle('') + guides(color=F) +
  scale_color_manual(values=plot_cols) +
  gtex_v8_figure_theme() +
  geom_hline(yintercept=0, linetype="dashed", colour="darkgrey") +
  theme(axis.text.x=element_text(angle=45,hjust=1),
        legend.title=element_blank(),
        legend.position=c(0.1,0.8),
        legend.key.height = unit(0.02, "cm"))

sfig <- plot_grid(sf1, sf2, nrow=2, labels=c('A', 'B'), axis='tlbr')

ggsave(sfig, file=paste0(data_dir, 'paper_figures/sfig_continuous_v2.pdf'), width=7.2, height=7, units="in")



