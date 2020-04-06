#!/usr/bin/env Rscript

## R script to pick tissues and individuals for imputation and outlier calling
## Run from upper level directory of the repo (or adjust path to tissue colors)

require(ggplot2)
require(gplots)
require(stringr)
require(data.table)

baseDir = Sys.getenv('RAREDIR')
dir = paste0(baseDir, '/preprocessing_v8')

##---------------- FUNCTIONS
## theme for functions and plots below
mytheme = theme_classic()

plot.ind.miss <- function(design, title = '', thresh = NULL) {
    ind.miss = data.frame(ind = colnames(design), perc = 1 - colMeans(design))
    ind.miss = ind.miss[order(ind.miss$perc), ]
    ind.miss.bar = ggplot(data = ind.miss, aes(x = c(1:nrow(ind.miss)), y = perc)) +
        geom_bar(stat = 'identity', fill = 'dodgerblue3')  +
        xlab('Individuals') + ylab('Missingness') + ggtitle(title) + mytheme
    if (!is.null(thresh)) {
        ind.miss.bar = ind.miss.bar + geom_hline(yintercept = thresh)
        keep = ind.miss$ind[ind.miss$perc <= thresh]
    } else {
        keep = ind.miss$ind
    }
    return(list(plot = ind.miss.bar, keep = keep))
}

plot.tissue.miss <- function(design, title = '', thresh = NULL) {
    tiss.miss = data.frame(tissue = rownames(design), perc = 1 - rowMeans(design))
    tiss.miss = tiss.miss[order(tiss.miss$perc), ]
    tiss.miss$tissue = factor(tiss.miss$tissue, levels = tiss.miss$tissue)
    tiss.miss.bar = ggplot(data = tiss.miss, aes(x = tissue, y = perc, fill = tissue)) +
        geom_bar(stat = 'identity') + coord_flip() + guides(fill = F) +
	xlab('') + ylab('Missingness') + ggtitle(title) +
        scale_fill_manual(values = gtex.colors) + mytheme
    if (!is.null(thresh)) {
        tiss.miss.bar = tiss.miss.bar + geom_hline(yintercept = thresh)
        keep = tiss.miss$tissue[tiss.miss$perc <= thresh]
    } else {
        keep = tiss.miss$tissue
    }
    return(list(plot = tiss.miss.bar, keep = keep))
}

##---------------- MAIN

## Read in official GTEx colors file
gtex.color.map = read.table('../gtex_tissue_colors.txt', sep = '\t', header = T, stringsAsFactors = F)
gtex.colors = paste0('#', gtex.color.map$tissue_color_hex)
names(gtex.colors) = gtex.color.map$tissue_site_detail_id

## Read in list of EA samples
eas.wgs = scan(paste0(dir, '/gtex_2017-06-05_v8_euro_VCFids.txt'), what = character())

## Read in sample to tissue correspondence and turn it into a individual to tissue correspondence
meta = read.table(paste0(dir, '/gtex_2017-06-05_v8_samples_tissues.txt'), header = F,
                  stringsAsFactors = F, col.names = c('Sample', 'Tissue'))
meta$Id = apply(str_split_fixed(meta$Sample, '-', 6)[, c(1:2)], 1, paste, collapse = '-')

## Read in normalized expression data and subset above correspondance to individuals with corrected data
## (this is a subset because it only includes individuals that were genotyped)
expr = read.table(gzfile(paste0(dir,'/gtex_2017-06-05_normalized_expression_v8ciseQTLs_removed.txt.gz')), header=T)
colnames(expr) = gsub("[.]", "-", colnames(expr))
meta = meta[meta$Id %in% colnames(expr), ]

## Get tissue and sample names
tissues = sort(unique(meta$Tissue))
brain.indices = grep('Brain', tissues)
individuals = unique(meta$Id)
gtex.colors = gtex.colors[tissues] # just to make sure they are ordered correctly

## Make matrix of samples by tissues with expression data
exp.design = matrix(0, ncol = length(individuals), nrow = length(tissues),
                    dimnames = list(tissues, individuals))
for (t in tissues) {
    inds = meta$Id[meta$Tissue == t]
    exp.design[t, inds] = 1
}

## Make heatmap of the samples by tissues
pdf(paste0(baseDir, '/figures/gtex_v8_design.pdf'), width = 10, height = 5)

heatmap.2(exp.design, Rowv = T, Colv = T, dendrogram = 'n', labCol = '', labRow = '',
          xlab = 'Individuals', ylab = 'Tissues', margins = c(2,2), lwid=c(0.1,4), lhei=c(0.7,4),
          RowSideColors = gtex.colors, col = c('white', 'dodgerblue3'), 
          trace = 'none', cexRow = .5, cexCol = .3, key = F, main = 'GTEx design matrix')

missingness.ind = plot.ind.miss(exp.design, 'Missingness by individual (original)', thresh = 0.75)
missingness.tiss = plot.tissue.miss(exp.design, 'Missingness by tissue (original)', thresh = 0.75)
#
print(missingness.ind$plot)
print(missingness.tiss$plot)
#
### Filter for tissues with at most 75% missingness
tissues.final = as.character(missingness.tiss$keep)
exp.design = exp.design[tissues.final, ]

## Filter for individuals with at most 75% missingness (based on subset set of tissues)
inds.final = as.character(plot.ind.miss(exp.design, thresh = 0.75)$keep)
exp.design = exp.design[, inds.final]

### Make heatmap and barplots with filtered tissues and individuals
col.ind = ifelse(colnames(exp.design) %in% eas.wgs, 'navyblue', 'white') # color column if individual has WGS
heatmap.2(exp.design, Rowv = T, Colv = T, dendrogram = 'n', labCol = '', labRow = '',
          xlab = 'Individuals', ylab = 'Tissues', margins = c(2,2), lwid=c(0.1,4), lhei=c(0.7,4),
          ColSideColors = col.ind,
          RowSideColors = gtex.colors[tissues.final], col = c('white', 'dodgerblue3'), 
          trace = 'none', cexRow = .5, cexCol = .3, key = F, main = 'GTEx design matrix, filtered')

print(plot.ind.miss(exp.design, 'Missingness by individual (filtered)')$plot)
print(plot.tissue.miss(exp.design, 'Missingness by tissue (filtered)')$plot)

dev.off()

## Calculate overall missingess in the filtered design matrix
cat('Average missingness', mean(exp.design), '\n')

## Write out selected tissues and individuals
write.table(sort(inds.final), paste0(dir, '/gtex_2017-06-05_v8_individuals_passed.txt'),
            sep = '\t', quote = F, col.names = F, row.names = F)
write.table(sort(tissues.final), paste0(dir, '/gtex_2017-06-05_v8_tissues_passed.txt'),
            sep = '\t', quote = F, col.names = F, row.names = F)
write.table(exp.design, paste0(dir, '/gtex_2017-06-05_v8_design_passed.txt'),
            sep = '\t', quote = F, col.names = T, row.names = T)

## Subset the normalized expression file to the individuals and tissues selected
expr.subset = expr[which(expr$Tissue %in% tissues.final),]
rm(expr)
expr.subset = expr.subset[, c('Tissue', 'Gene', sort(inds.final))]

## Then subset to the genes expressed in each tissue that are either protein coding or lincRNA
## first reading in the list of autosomal protein-coding  & lincRNA genes
autosomal.df = read.table(paste0(dir, '/gencode.v26.GRCh38.genes_genetypes_autosomal_PCandlinc_only.txt'),
                          header = F, stringsAsFactors = F)
autosomal.selected = autosomal.df[autosomal.df[, 2] %in% c('lincRNA', 'protein_coding'), 1]

print(length(unique(expr.subset$Gene)))
genes.counts = table(expr.subset$Gene)
genes.keep = names(genes.counts)[which(genes.counts == length(tissues.final))]
genes.keep = genes.keep[which(genes.keep %in% autosomal.selected)]
print(length(genes.keep))
expr.subset = expr.subset[which(expr.subset$Gene %in% genes.keep), ]

## finally restandardize and output subsetted expression matrix
expr.subset = as.data.frame(expr.subset)
expr.subset[, 3:ncol(expr.subset)] = t(scale(t(expr.subset[, 3:ncol(expr.subset)])))

gzout = gzfile(paste0(dir, '/gtex_2017-06-05_normalized_expression_v8ciseQTLs_removed.subset.txt.gz'), 'w')
write.table(expr.subset, gzout, sep = '\t', quote = F, col.names = T, row.names = F)
close(gzout)
