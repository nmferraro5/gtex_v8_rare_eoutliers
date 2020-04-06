library(data.table)
library(dplyr)
library(ggplot2)
library(FNN)

data_dir = '/Users/nicoleferraro/durga_local/data/goats_data/v8_data/mpra/'

coloc_test = fread(paste0(data_dir, 'mpra_coloc_test_variants.txt'))
coloc_test$af_gnomad = sapply(coloc_test$af_gnomad, function(x) ifelse(is.na(x), 0, x))
coloc_control = fread(paste0(data_dir, 'mpra_coloc_test_variants_all_controls.txt'))
coloc_control$af_gnomad = sapply(coloc_control$af_gnomad, function(x) ifelse(is.na(x), 0, x))

i=0
notSame = TRUE
while(notSame) {
  i = i + 1
  print(i)
  sampled_controls = coloc_control %>% sample_n(50)
  af_pval = wilcox.test(sampled_controls$af_gnomad, coloc_test$af_gnomad)$p.value
  dist_pval = wilcox.test(sampled_controls$distTSS, coloc_test$distTSS)$p.value
  if (af_pval > 0.1 && dist_pval > 0.1) {
    notSame = FALSE
  } else if (i > 1000) {
    notSame = FALSE
  }
}

both_values = rbind(sampled_controls %>% select(distTSS, af_gnomad) %>% mutate(Cat = 'Control'), 
                    coloc_test %>% select(distTSS, af_gnomad) %>% mutate(Cat = 'Test')) 

ggplot(both_values, aes(x=Cat, y=af_gnomad)) + geom_boxplot() + theme_bw() + xlab('')
ggplot(both_values, aes(x=Cat, y=distTSS)) + geom_boxplot() + theme_bw() + xlab('')

write.table(sampled_controls, file=paste0(data_dir, 'mpra_coloc_test_variants_sampled_controls.txt'), sep='\t', quote=F, row.names=F)

ws_test = fread(paste0(data_dir, 'mpra_watershed_test_variants_v2.txt'))
ws_test$af_gnomad = sapply(ws_test$af_gnomad, function(x) ifelse(is.na(x), 0, x))
ws_test = ws_test %>% filter(!(variant_cat %in% c('coding', 'splice', 'stop'))) %>%
  select(Gene,Chr,Pos,Ref,Allele,af_gnomad,distTSS) %>% distinct()
ws_control = fread(paste0(data_dir, 'mpra_watershed_test_variants_all_controls_v2.txt'))
ws_control$af_gnomad = sapply(ws_control$af_gnomad, function(x) ifelse(is.na(x), 0, x))
ws_control = ws_control %>% filter(!(variant_cat %in% c('coding', 'splice', 'stop'))) %>%
  select(Gene,Chr,Pos,Ref,Allele,af_gnomad,distTSS) %>% distinct()

get_match <- function(i) {
  tdata = ws_test[i,]
  tcontrols = filter(ws_control, Gene == tdata$Gene[1])
  if (nrow(tcontrols) == 0) {
    return(rep(NA, ncol(tcontrols)))
  }
  if (nrow(tcontrols) == 1) {
    tneighbor = get.knnx(tcontrols %>% select(af_gnomad,distTSS), tdata %>% select(af_gnomad,distTSS), k=1)
  } else {
    tneighbor = get.knnx(tcontrols %>% select(af_gnomad,distTSS), tdata %>% select(af_gnomad,distTSS), k=2)
  }
  tcontrol = tcontrols[tneighbor$nn.index,]
  return(tcontrol)
}

chosen_controls = do.call(rbind, lapply(1:nrow(ws_test), function(x) get_match(x))) %>%
  filter(!is.na(Gene))
#set.seed(92413) # used 92413 for v1, and chose distance based on af_gnomad
set.seed(1412) # used 92413 for v1
ws_sampled_controls = chosen_controls %>% distinct() %>% sample_n(1000)
wilcox.test(ws_sampled_controls$af_gnomad, ws_test$af_gnomad)
wilcox.test(ws_sampled_controls$distTSS, ws_test$distTSS)

both_values = rbind(ws_sampled_controls %>% select(distTSS, af_gnomad) %>% mutate(Cat = 'Control'), 
                    ws_test %>% select(distTSS, af_gnomad) %>% mutate(Cat = 'Test')) 

ggplot(both_values, aes(x=Cat, y=af_gnomad)) + geom_boxplot() + theme_bw() + xlab('')
ggplot(both_values, aes(x=Cat, y=distTSS)) + geom_boxplot() + theme_bw() + xlab('')

write.table(sampled_controls, file=paste0(data_dir, 'mpra_watershed_test_variants_sampled_controls.txt'), sep='\t', quote=F, row.names=F)

test_moderate = fread(paste0(data_dir, 'mpra_watershed_test_variants_moderate.txt'))
test_moderate$af_gnomad = sapply(test_moderate$af_gnomad, function(x) ifelse(is.na(x), 0, x))
test_moderate = test_moderate %>% filter(!(variant_cat %in% c('coding', 'splice', 'stop'))) %>%
  select(Gene,Chr,Pos,Ref,Allele,af_gnomad,distTSS) %>% distinct()

get_match <- function(i) {
  tdata = ws_test[i,]
  tcontrols = filter(test_moderate, Gene == tdata$Gene[1])
  if (nrow(tcontrols) == 0) {
    return(rep(NA, ncol(tcontrols)))
  }
  if (nrow(tcontrols) < 2) {
    tneighbor = get.knnx(tcontrols %>% select(af_gnomad), tdata %>% select(af_gnomad), k=1)
  } else {
    tneighbor = get.knnx(tcontrols %>% select(af_gnomad), tdata %>% select(af_gnomad), k=2)
  }
  tcontrol = tcontrols[tneighbor$nn.index,]
  return(tcontrol)
}

test_moderate = filter(test_moderate, distTSS < max(ws_test$distTSS))
test_moderate = filter(test_moderate, distTSS > min(ws_test$distTSS))
chosen_moderates = do.call(rbind, lapply(1:nrow(ws_test), function(x) get_match(x))) %>%
  filter(!is.na(Gene))

kseeds=sample(1:50000000, 2000)
good_seeds = c()
for (ks in kseeds) {
  set.seed(ks) #521 was good
  sampled_moderate = chosen_moderates %>% distinct() %>% sample_n(1024)
  p1 = wilcox.test(sampled_moderate$af_gnomad, ws_test$af_gnomad)$p.value
  p2 = wilcox.test(sampled_moderate$distTSS, ws_test$distTSS)$p.value
  if (p1 > 0.05 & p2 > 0.05) {
    good_seeds = c(good_seeds, ks)
  }
}
good_seeds # use 27499160

both_values = rbind(sampled_moderate %>% select(distTSS, af_gnomad) %>% mutate(Cat = 'Moderate'), 
                    ws_test %>% select(distTSS, af_gnomad) %>% mutate(Cat = 'Test')) 

ggplot(both_values, aes(x=Cat, y=af_gnomad)) + geom_boxplot() + theme_bw() + xlab('')
ggplot(both_values, aes(x=Cat, y=distTSS)) + geom_boxplot() + theme_bw() + xlab('')

write.table(sampled_controls, file=paste0(data_dir, 'mpra_watershed_test_variants_sampled_controls.txt'), sep='\t', quote=F, row.names=F)



