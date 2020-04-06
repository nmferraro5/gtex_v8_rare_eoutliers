#!/usr/bin/env Rscript

## From the multi-tissue outliers, determine which tissues have extreme effects

library(data.table)
library(stringr)
library(optparse)
library(ggplot2)
library(dplyr)
library(viridis)

dir = Sys.getenv('RAREDIR')


source('../paper_figures/plot_settings.R')

##################### FUNCTIONS

## Function that returns the tissues or data type with |Z-score| > 2
## Assumes that a data.table called expr with the expression data as constructed in preprocessing/ exists and is available
## Input: vector with Ind and Gene that defines the outlier, in columns 1 and 2, respectively
##        threshold for calling a tissue or data type extreme
##        name of the column with the tissue or data type
## Output: comma-separated string of tissues or data type with extreme expression
get.extreme <- function(outlier.info, thresh = 3, col.name = 'Tissue', mz=FALSE, restrict=FALSE) {
    ind = outlier.info[1]
    gene = outlier.info[2]
    expr.df = as.data.frame(expr[gene, c(col.name, ind), with = F])
    #expr.df = as.data.frame(expr[gene, c(col.name, ind)])
    if (mz == TRUE) {
      rind = which(abs(expr.df[,2]) == max(abs(expr.df[,2]),na.rm=T))
      #rind = which(expr.df[,2] == min(expr.df[,2],na.rm=T))
      return(expr.df[rind,2])
    }
    indices = abs(expr.df[, ind]) > thresh & !is.na(expr.df[, ind])
    #indices = expr.df[, ind] < 2*pnorm(-abs(thresh)) & !is.na(expr.df[, ind])
    # Adding these 2 lines to restrict to tissue-specific if remainder are < 1
    if (length(which(indices)) > 0 & restrict) {
        indices = (abs(expr.df[, ind]) > thresh & !is.na(expr.df[, ind])) | (abs(expr.df[, ind]) < thresh & abs(expr.df[, ind]) > 2 & !is.na(expr.df[, ind]))
    }
    tissues = sort(expr.df[indices, col.name]) # sorted so order will match ER map file (for GTEx data)
    return(paste(tissues, collapse = ','))
}

## Function to compute the jaccard index for pairs of tissues among the list of extreme tissues.
## Input: the name of the outlier calling method, outlier list
## Output: Data frame with values from an upper triangular matrix and associated labels and positions.
compute.jaccard.pairs <- function(method, outliers.list) {
    extreme.method = outliers.list %>% filter(Method == method) %>% select(Tissues) %>% unlist
    ntis = length(tissues.ordered)
    nrows = ntis * (ntis-1) / 2
    norm.df = data.frame(t1 = character(nrows), t2 = character(nrows),
                         x = numeric(nrows), y = numeric(nrows),
                         val = numeric(nrows), stringsAsFactors = FALSE)
    row = 1
    for (i in 1:(ntis-1)) {
        for (j in (i+1):ntis) {
            norm.df[row, ] = list(tissues.ordered[i], tissues.ordered[j], i, j,
                                  jaccard(extreme.method, tissues.ordered[i], tissues.ordered[j]))
            row = row + 1
        }
    }
    norm.df = norm.df %>% mutate(t1 = factor(t1, levels = tissues.ordered),
                                 t2 = factor(t2, levels = tissues.ordered))
    return(norm.df)
}

## Function to compute the Jaccard index for two given tissues.
## Input: The vector of extreme tissues, the names of the two tissues.
## The Jaccard index (intersection/union) of the occurence of the tissues among the extreme tissues list
jaccard <- function(extreme.tissues, tissue1, tissue2) {
    indices1 = grep(tissue1, extreme.tissues, fixed = TRUE)
    indices2 = grep(tissue2, extreme.tissues, fixed = TRUE)
    jaccard = length(intersect(indices1, indices2)) / length(union(indices1, indices2))
    return(jaccard)
}

## Function to generate point plot of jaccard distances
## Input: data frame returned by compute jaccard pairs
## Note: labels are cut off. This can be fixed in the future
plot.pairs <- function(pairs, method) {
    p = ggplot(pairs, aes(x = t1, y = t2)) +
        geom_point(aes(size = val, colour = val)) +
        geom_text(data = filter(pairs, x == y-1), aes(label = t1),
                  nudge_x = -0.1, nudge_y = -0.4, hjust = 0, size = 3, angle = -30) +
        xlab('') + ylab('') +
        scale_color_viridis(name = 'Jaccard') +
        guides(size = FALSE) +
        theme_bw() + ggtitle(method) +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank())
    return(p)
}

##################### MAIN

option_list = list(make_option(c('--Z.SCORES'), type = 'character', default = NULL, help = 'path to the Z-score data'),
        make_option(c('--EXP.DATA'), type = 'character', default = NULL, help = 'path to gzipped expression data file'))

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

zscore_file = opt$Z.SCORES
exp_file = opt$EXP.DATA

## Get outliers from top outlier/control data frames and combine into a single data frame
medz_outliers = fread(zscore_file) %>% mutate(Method = 'MEDZ')
outlier_list_top = list(medz_outliers)
outliers.top = do.call(rbind, lapply(outlier_list_top, function(df) df[df$Y == 'outlier', c('Ind','Gene','Method')]))

## Load expression data
expr = fread(paste0('zcat ', exp_file))
setkey(expr, Gene)

outliers.top$Tissues = apply(outliers.top, 1, get.extreme)
outliers.top$Nextreme = str_count(outliers.top$Tissues, ',') + 1 # +1 because fence post problem

## For outliers with extreme expression in 3 or fewer tissues, are the tissues related ?
outliers.specific = filter(outliers.top, Nextreme <= 3) %>%
    mutate(Method = factor(Method, levels = names(meth.colors)))

plot.spec = ggplot(outliers.specific, aes(x = Nextreme, fill = Method)) +
    geom_bar(position = position_dodge(width = 0.9)) +
    xlab('Number of tissues with |Z-score| > 3') +
    scale_fill_manual(values = meth.colors) + mytheme

## Reorganize data for hierarchical clustering then run clustering
genes = unique(expr[, Gene])
tissues = unique(expr[, Tissue])

expr.reshaped = matrix(NA, nrow = length(tissues), ncol = (ncol(expr) - 2) * length(genes))
rownames(expr.reshaped) = tissues
for (t in tissues) {
    expr.reshaped[t,] = unlist(expr[Tissue == t, c(3:ncol(expr))])
}

tissue.clustering = hclust(dist(expr.reshaped), method = 'average')

tissues.ordered = tissues[tissue.clustering$order]

## for each method, get a matrix of the coocurrence of tissue pairs among outliers
## plot this
pdf(paste0(dir, '/figures/GTEXv8_pair_jaccard.pdf'), width = 10, height = 8)

print(plot.spec)

plot(tissue.clustering)

for (meth in unique(outliers.top$Method)) {
    pair.values = compute.jaccard.pairs(meth, outliers.top)
    print(plot.pairs(pair.values, meth))
}

for (meth in unique(outliers.top$Method)) {
    pair.values.spec = compute.jaccard.pairs(meth, filter(outliers.top, Nextreme == 2))
    print(plot.pairs(pair.values.spec, paste(meth, 'tissue specific')))
}

dev.off()
