library(data.table)
library(dplyr)

data_dir = '/srv/scratch/restricted/GOATs/data_v8/splicing/sqtl/formatted/'

enloc_files = dir(data_dir, '*.enloc.rst.gz')

process_file <- function(ef) {
  print(ef)
  edata = fread(paste0('zcat ', data_dir, ef),data.table=F,fill=T)
  edata = filter(edata, locus_rcp > 0.5)
  if (nrow(edata) > 0) {
    trait = strsplit(ef, '__PM__')[[1]][1]
    tissueRest = strsplit(ef, '__PM__')[[1]][2]
    tissue = strsplit(tissueRest, '.enloc.rst.gz')[[1]][1]
    edata$Trait = trait
    edata$Tissue = tissue
    return(edata)
  }
}

all_enloc_data = do.call(rbind, lapply(enloc_files, process_file))

write.table(all_enloc_data, file=paste0(data_dir, 'sqtl_enloc_rcp0.5_collapsed_update.txt'),sep='\t',quote=F,row.names=F)

