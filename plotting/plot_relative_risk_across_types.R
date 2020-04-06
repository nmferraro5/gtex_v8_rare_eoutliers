#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
require(dplyr)
library(reshape)
library(scales)
require(RColorBrewer)
library(epitools)

#data_dir = '/Users/nicoleferraro/durga_local/data/goats_data/v8_data/'
data_dir = '/users/nferraro/data/goats_data/v8_data/'

### splicing
print('Reading splicing')
splice_data = fread(paste0('zcat ', data_dir, 'splicing/gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.splice.variantAnnotations.filtered.txt.gz'),data.table=F)
splice_data = splice_data %>% select(indiv_id,gene_id,value,pval_bin,tier2,af_gtex,categoryOutlier,variant_cat,sv_v7,af_gnomad,af_bin,color_R,color_G,color_B)
splice_data = filter(splice_data,sv_v7==1)
new_cats = sapply(1:nrow(splice_data), function(x) ifelse(splice_data$variant_cat[x] == 'splice', splice_data$tier2[x], splice_data$variant_cat[x]))
splice_data$variant_cat = new_cats
splice_data$OutlierValue = -log10(splice_data$value)
mval = max(filter(splice_data, is.finite(OutlierValue))$OutlierValue)
splice_data$OutlierValue = sapply(splice_data$OutlierValue, function(x) ifelse(is.finite(x), x, mval))

### expression
print('Reading expression')
exp_data = fread(paste0('zcat ', data_dir, 'gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.medz.variantAnnotations.filtered.txt.gz'),data.table=F)
exp_data = exp_data %>% select(indiv_id,gene_id,MedZ,Y,tier2,af_gtex,categoryOutlier,variant_cat,sv_v7,af_gnomad,medz_bin,color_R,color_G,color_B)
exp_data = filter(exp_data,sv_v7==1)
new_cats = sapply(1:nrow(exp_data), function(x) ifelse(exp_data$variant_cat[x] == 'splice', exp_data$tier2[x], exp_data$variant_cat[x]))
exp_data$variant_cat = new_cats
exp_data$OutlierValue = -log10(2*pnorm(-abs(exp_data$MedZ)))

### ASE
print('Reading ASE')
ase_data = fread(paste0('zcat ', data_dir, 'ASE/gtexV8.outlier.controls.v8ciseQTLs.globalOutliers.removed.ase.variantAnnotations.v8.vg.txt.gz'),data.table=F)
ase_data = ase_data %>% select(indiv_id,gene_id,value,pval_bin,tier2,variant_cat,sv_v7,af_gnomad,variant_color1,variant_color2,variant_color3)
ase_data = filter(ase_data,sv_v7==1)
new_cats = sapply(1:nrow(ase_data), function(x) ifelse(ase_data$variant_cat[x] == 'splice', ase_data$tier2[x], ase_data$variant_cat[x]))
ase_data$variant_cat = new_cats
ase_data$OutlierValue = -log10(ase_data$value)

### get relative risk per category
splice_outliers = filter(splice_data, value < 0.0027)
splice_controls = filter(splice_data, value > 0.0027)
exp_outliers = filter(exp_data, abs(MedZ) > 3)
exp_controls = filter(exp_data, abs(MedZ) < 3)
ase_outliers = filter(ase_data, value < 0.0027)
ase_controls = filter(ase_data, value > 0.0027)

risks = data.frame(Risk = numeric(), Lower = numeric(), Upper = numeric(), Pval = numeric(), St = character())
vcats = unique(exp_data$variant_cat)
### Relative risk
# for (vcat in vcats) {
#   print(vcat)
#   splice_nn = nrow(filter(splice_controls, variant_cat != vcat))
#   splice_ny = nrow(filter(splice_controls, variant_cat == vcat))
#   splice_yn = nrow(filter(splice_outliers, variant_cat != vcat))
#   splice_yy = nrow(filter(splice_outliers, variant_cat == vcat))
#   splicetable = rbind(c(splice_nn,splice_ny),c(splice_yn,splice_yy))
#   srr = epitab(splicetable, method = 'riskratio')
#   risks = rbind(risks, data.frame(Risk = srr$tab[2,5],
#                                   Lower = srr$tab[2,6],
#                                   Upper = srr$tab[2,7],
#                                   Pval = srr$tab[2,8],
#                                   Cat = vcat,
#                                   Type = 'Splicing'))
# 
#   exp_nn = nrow(filter(exp_controls, variant_cat != vcat))
#   exp_ny = nrow(filter(exp_controls, variant_cat == vcat))
#   exp_yn = nrow(filter(exp_outliers, variant_cat != vcat))
#   exp_yy = nrow(filter(exp_outliers, variant_cat == vcat))
#   exptable = rbind(c(exp_nn,exp_ny),c(exp_yn,exp_yy))
#   err = epitab(exptable, method = 'riskratio')
#   risks = rbind(risks, data.frame(Risk = err$tab[2,5],
#                                   Lower = err$tab[2,6],
#                                   Upper = err$tab[2,7],
#                                   Pval = err$tab[2,8],
#                                   Cat = vcat,
#                                   Type = 'Total expression'))
#   
#   ase_nn = nrow(filter(ase_controls, variant_cat != vcat))
#   ase_ny = nrow(filter(ase_controls, variant_cat == vcat))
#   ase_yn = nrow(filter(ase_outliers, variant_cat != vcat))
#   ase_yy = nrow(filter(ase_outliers, variant_cat == vcat))
#   asetable = rbind(c(ase_nn,ase_ny),c(ase_yn,ase_yy))
#   arr = epitab(asetable, method = 'riskratio')
#   risks = rbind(risks, data.frame(Risk = arr$tab[2,5],
#                                   Lower = arr$tab[2,6],
#                                   Upper = arr$tab[2,7],
#                                   Pval = arr$tab[2,8],
#                                   Cat = vcat,
#                                   Type = 'ASE'))
# }

all_coefs = data.frame(Beta = numeric(), SE = numeric(), Zval = numeric(), Pval = numeric(), Category = character(), Variant = character())

### Continuous risk
for (vcat in vcats) {
  print(vcat)
  splice_data = splice_data %>% mutate(HasCAT = ifelse(variant_cat == vcat, 1, 0))
  splice_lm = as.data.frame(summary(glm(HasCAT ~ OutlierValue, data = splice_data, family = binomial))$coefficients)[2,]
  colnames(splice_lm) = c('Beta', 'SE', 'Zval', 'Pval')

  ase_data = ase_data %>% mutate(HasCAT = ifelse(variant_cat == vcat, 1, 0))
  ase_lm = as.data.frame(summary(glm(HasCAT ~ OutlierValue, data = ase_data, family = binomial))$coefficients)[2,]
  colnames(ase_lm) = c('Beta', 'SE', 'Zval', 'Pval')
  
  exp_data = exp_data %>% mutate(HasCAT = ifelse(variant_cat == vcat, 1, 0))
  exp_lm = as.data.frame(summary(glm(HasCAT ~ OutlierValue, data = exp_data, family = binomial))$coefficients)[2,]
  colnames(exp_lm) = c('Beta', 'SE', 'Zval', 'Pval')
  
  new_coefs = rbind(splice_lm %>% mutate(Category = 'sOutliers', Variant = vcat),
                    exp_lm %>% mutate(Category = 'eOutliers', Variant = vcat),
                    ase_lm %>% mutate(Category = 'aseOutliers', Variant = vcat))
  all_coefs = rbind(all_coefs, new_coefs)
}

# risks = risks %>% arrange(by=Risk) 
# risks$Cat = factor(risks$Cat, levels=unique(risks$Cat))
# risks$Type = factor(risks$Type, levels=c('ASE','Splicing', 'Total expression'))

pcols = exp_data %>% filter(!(medz_bin %in% c(NA,"0.000~0.002"))) %>%
  group_by(variant_cat) %>% sample_n(1) %>%
  ungroup() %>% mutate(CatCol = rgb(color_R, color_G, color_B)) %>%
  select(variant_cat, CatCol)
plot_cols = c(pcols$CatCol[1:11],'#3e6690','#33669a',pcols$CatCol[14:16])
names(plot_cols) = unique(pcols$variant_cat)

#save(risks, plot_cols, file=paste0(data_dir, 'gtexV8.relative.risks.all.types.RData'))
save(all_coefs, plot_cols, file=paste0(data_dir, 'gtexV8.continuous.risks.all.types.RData'))


load(paste0(data_dir, 'gtexV8.relative.risks.all.types.RData'))
risks$Cat = factor(risks$Cat, levels=c('no_variant','other_noncoding','TE', 'coding','TSS', 'conserved_noncoding','INV','BND','DEL','CNV','DUP', 'splice_region_variant', 'splice_acceptor_variant','frameshift','splice_donor_variant', 'stop'))
plot_cols['no_variant'] = 'darkgrey'

ggplot(risks %>% filter(Risk < 70), aes(x=Cat,y=Risk,Group=Type)) +
  geom_point(size=5, aes(color=Cat,shape=Type),position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width=0,position=position_dodge(width=0.5)) +
  theme_bw() + ylab('') + xlab('') + scale_y_continuous(trans='log2') +
  ggtitle('') + guides(color=F,shape=F) +
  scale_color_manual(values=plot_cols) +
  geom_hline(yintercept=1, linetype="dashed", colour="darkgrey") +
  theme(axis.title = element_text(size=20),
        axis.text.y = element_text(size=16),
        axis.text.x=element_text(size=14,angle=45,hjust=1),
        legend.title=element_blank(),
        legend.text=element_text(size=20)) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

ggplot(risks %>% filter(Risk < 70), aes(x=Cat,y=Risk,Group=Type)) +
  geom_point(size=5, aes(color=Cat,shape=Type),position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width=0,position=position_dodge(width=0.5)) +
  theme_bw() + ylab('') + xlab('') + scale_y_continuous(trans='log2') +
  ggtitle('') + guides(color=F,shape=F) +
  scale_color_manual(values=plot_cols) +
  geom_hline(yintercept=1, linetype="dashed", colour="darkgrey") +
  theme(axis.title = element_text(size=20),
        axis.text.y = element_text(size=22),
        axis.text.x=element_text(size=18,angle=45,hjust=1),
        legend.title=element_blank(),
        legend.text=element_text(size=20)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

load(paste0(data_dir, 'gtexV8.continuous.risks.all.types.RData'))
all_coefs$Variant = factor(all_coefs$Variant, levels=c('no_variant','other_noncoding','TE', 'coding','TSS', 'conserved_noncoding','INV','BND','DEL','CNV','DUP', 'splice_region_variant', 'splice_acceptor_variant','frameshift','splice_donor_variant', 'stop'))
plot_cols['no_variant'] = 'darkgrey'

ggplot(all_coefs, aes(x=Variant,y=Beta,Group=Category)) +
  geom_point(size=5, aes(color=Variant,shape=Category),position=position_dodge(width=0.5)) +
  geom_errorbar(aes(ymin = Beta - 1.95*SE, ymax = Beta + 1.95*SE), width=0,position=position_dodge(width=0.5)) +
  theme_bw() + ylab('') + xlab('') + ylim(c(-0.5,1.5)) +
  ggtitle('') + guides(color=F,shape=F) +
  scale_color_manual(values=plot_cols) +
  geom_hline(yintercept=0, linetype="dashed", colour="darkgrey") +
  theme(axis.title = element_text(size=20),
        axis.text.y = element_text(size=16),
        axis.text.x=element_text(size=14,angle=45,hjust=1),
        legend.title=element_blank(),
        legend.text=element_text(size=20)) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))



