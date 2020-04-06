#!/usr/bin/env Rscript

## Helper functions for compute_enrichments.R

## Load required packages
library(dplyr)
library(data.table)
library(epitools)

#### Function to select outliers based on a significance or effect size threshold
## Input: File of outlier test statistics (assumes the file name is period-separated and ends in <method>.txt), 
## the significance threshold for defining outliers (or number of outliers to keep)
## if topn is true, then keep the top threshold number of outliers, otherwise use threshold as a significance threshold
## Output: Assignments of outliers and controls for genes with outliers
filter.outliers <- function(stats_file, threshold, topn = FALSE) {
    print(stats_file)
    test.stats = fread(stats_file)
    outliers = test.stats %>% filter(Y == 'outlier')
    # filename must end in <method>.txt where method, is the name of the method in question
    method = toupper(gsub('.*\\.','', sub('.txt', '', basename(stats_file))))
    if (method == 'MEDZ') {
        ## filter by medz
        outliers = outliers %>% mutate(Rank = -abs(MedZ)) # make negative, so order is from smallest to biggest
        threshold = -abs(threshold) # so that thresholding will work on transformed ranking, works for topn too - see below
    } else {
        ## filter by FDR
        outliers = outliers %>% mutate(Rank = FDR)
    }

    ## use ranking to subset outliers
    if (topn) {
        outliers = outliers %>% top_n(-abs(threshold), Rank) %>% mutate(Y='outlier') %>% mutate(UID=paste0(Ind,Gene)) # negative first argument to get smallest n
        controls = filter(test.stats, Gene %in% outliers$Gene) %>% mutate(UID=paste0(Ind,Gene)) %>% filter(!(UID %in% outliers$UID)) %>% mutate(Y='control')
        outliers = outliers %>% select(-UID,-Rank)
        controls = controls %>% select(-UID)
    } else {
        outliers = outliers %>% filter(Rank < threshold) %>% mutate(Y='outlier')
        controls = filter(test.stats,Gene %in% outliers$Gene,abs(MedZ)<= abs(threshold)) %>% mutate(Y='control')
        outliers = outliers %>% select(-Rank)
    }
    
    ## Only return genes with at least one outlier called
    #test.stats = test.stats %>% filter(Gene %in% outliers$Gene) %>% mutate(Method = method)
    test.stats = rbind(outliers,controls) %>% mutate(Method=method)
    return(test.stats)
}

## Merge features for a given individual for each gene in outlier/control data frame
## Note that everything after the period gets stripped off the gene ids for compatibility with datasets that exclude that part
## Input: individual ID, path to feature files, suffix of feature files, outlier/control data frame
## Output: Merged data frame including outlier info and features
get.features <- function(individual, feature_dir, feature_suffix, outliers) {
    outliers = outliers %>%
        filter(Ind == individual) %>%
        mutate(Y = ifelse(Y == 'outlier', 1, 0), # codify outlier status in preparation for logisitic regression
               Gene = sub('\\.[0-9]+','', Gene)) # strip part of gene id after period
    if ('Specific.Group' %in% colnames(outliers)) {
        outliers = select(outliers, Gene, Ind, Y, value, Specific.Group, Maf, Type) # removed Method, FDR, and medz_value
      
    } else {
        outliers = select(outliers, Gene, Ind, Y, Method, Maf, Type)
    }
    genes = outliers %>% select(Gene) %>% unlist
    indiv_file = paste0(feature_dir, individual, feature_suffix)
    features = read.table(indiv_file, sep = '\t', header = TRUE, fill=T) %>%
        mutate(Gene = sub('\\.[0-9]+','', gene_id)) %>% # same as above
        select(-gene_id) %>%
        filter(Gene %in% genes)
    ## merge outliers and features based on gene (both are already subset to the correct individual)
    outliers_with_fts = merge(outliers, features, by = 'Gene')
    return(outliers_with_fts)
}

### Get data frame with merged outliers and features across all individuals.
### Wrapper for get.features
### Note: outliers is listed as the first argument so it is easy to lapply across outlier sets
get.all.features <- function(outliers, maf, vartype, prefix, suffix) {
    featdir = paste0(prefix, maf, '/', vartype, '/')
    outliers = outliers %>% mutate(Maf = maf, Type = vartype)
    filenames = dir(featdir, suffix)
    inds_feats = sub(suffix, '', filenames, fixed = T)
    overlap_inds = intersect(inds_feats, outliers$Ind)
    outliers_feats = do.call(rbind, lapply(overlap_inds, get.features, featdir, suffix, outliers))
    return(outliers_feats)
}

## Get enrichments for all features for the given maf and variant type
## Input: MAF, variant type, outlier/control data frame, prefix of path to features,
##   suffix of feature file, whether or not to scale features, wherther they are count features
## Output: Data frame with estimates from logistic regreassion for each feature
compute.enrichments <- function(maf, vartype, outliers, prefix, suffix, scale, counts = FALSE) {
    ## Get merged outliers with features
    method = outliers$Method[1]
    cat(method, maf, vartype, '\n')
    outliers_with_feats = get.all.features(outliers, maf, vartype, prefix, suffix)
    ## if counts features, then create and 'any' feature and change the name of the exisiting feature
    ## (so it doesn't clash with the n_variants of the full feature set)
    if (counts) {
        outliers_with_feats = outliers_with_feats %>%
            mutate(n_variants_all = n_variants,
                   any_variant = as.numeric(n_variants_all > 0)) %>%
            select(-n_variants)
    } else {
        cadd.indices = colnames(outliers_with_feats)[grep('CADD', colnames(outliers_with_feats), fixed = T)]
        for (ind in cadd.indices) {
            outliers_with_feats[, ind] = ifelse(is.na(outliers_with_feats[, ind]), 0, outliers_with_feats[, ind])
        }
    }
    ## Call logistic regressions for each of the features and gather log odds
    if ('Specific.Group' %in% colnames(outliers_with_feats)) {
        start_ind = 8
    } else {
        start_ind = 7
    }
    features = colnames(outliers_with_feats)[start_ind:length(colnames(outliers_with_feats))]        
    estims = data.frame(Beta = numeric(), Stderr = numeric(), Pval = numeric(), Feature = character())
    for (feature in features) {
        ## Deal with case where all values for a given feature may be NA, or the variance is too small
        featurevar = var(outliers_with_feats[, feature], na.rm = TRUE)
        if (!is.na(featurevar) && featurevar > 0) {
            if (scale) {
                lmod = summary(glm(Y ~ scale(outliers_with_feats[, feature]), data = outliers_with_feats, family = binomial))
            } else {
                lmod = summary(glm(Y ~ outliers_with_feats[, feature], data = outliers_with_feats, family = binomial))
            }
            estims = rbind(estims, data.frame(Beta = lmod$coefficients[2, 1],
                                              Stderr = lmod$coefficients[2, 2],
                                              Pval = lmod$coefficients[2, 4],
                                              Feature = feature))
        } else {
            cat('MAF:', maf, 'vartype:', vartype, 'Method:', method, 'Feature:', feature, 'Feature var:', featurevar, '\n')
        }
    }
    ## Include the method and variant type in return data frame for ease of plotting
    estims = estims %>% mutate(Method = method, Type = vartype, Maf = maf)
    return(estims)
}

compute.relative.risk <- function(maf, vartype, outliers, prefix, suffix, scale, counts = TRUE) {
  ## Get merged outliers with features
  method = outliers$Method[1]
  cat(method, maf, vartype, '\n')
  if (vartype == 'HallLabSV') {
    outliers_with_feats = get.all.features(outliers, maf, vartype, prefix, suffix)
  } else {
    #outliers_with_feats = get.all.features(outliers, maf, vartype, prefix, suffix)
    vartype='SNPs+indels' # added to combine SNPs and indels
    outliers_with_feats = do.call(rbind, lapply(c('SNPs','indels'), function(x)
      get.all.features(outliers, maf, x, prefix, suffix)))
    outliers_with_feats = outliers_with_feats %>% group_by(Gene,Ind) %>%
      top_n(1,n_variants) %>% sample_n(1) %>% ungroup()
  }
  ## if counts features, then create and 'any' feature and change the name of the exisiting feature
  ## (so it doesn't clash with the n_variants of the full feature set)
  ## if the full feature set, zero-impute NAs for CADD TFBS features
  outliers_with_feats = outliers_with_feats %>%
    mutate(n_variants_all = n_variants,
           any_variant = as.numeric(n_variants_all > 0)) %>%
    select(-n_variants,-n_variants_all)
  ## Call logistic regressions for each of the features and gather log odds
  if ('Specific.Group' %in% colnames(outliers_with_feats)) {
    start_ind = 8
  } else {
    start_ind = 7
  }
  features = colnames(outliers_with_feats)[start_ind:length(colnames(outliers_with_feats))]
  estims = data.frame(Riskratio = numeric(), Lower = numeric(), Upper = numeric(), Pval = numeric(), Feature = character())
  for (feature in features) {
    ## Deal with case where all values for a given feature may be NA, or the variance is too small
    featurevar = var(outliers_with_feats[, feature], na.rm = TRUE)
    if (!is.na(featurevar) && featurevar > 0) {
      #lmod = summary(glm(Y ~ outliers_with_feats[, feature], data = outliers_with_feats, family = binomial))
      outliers = filter(outliers_with_feats, Y == 1)
      controls = filter(outliers_with_feats, Y == 0)
      outliers$Y = factor(outliers$Y, levels=c(0,1))
      controls$Y = factor(controls$Y, levels=c(0,1))
      outliers = as.data.frame(outliers)
      controls = as.data.frame(controls)
      print(length(unique(outliers[,feature])))
      if (length(unique(outliers[,feature])) > 1) {
        counttable = rbind(table(controls[,feature]),
                           table(outliers[,feature]))
      } else {
        counttable = rbind(table(controls[,feature]),
                           c(table(outliers[,feature]),0))
      }
      tryCatch({
        rr = epitab(counttable,method="riskratio")$tab
      }, error = function(e) {
        print(table(controls[,feature]))
        print(table(outliers[,feature]))
        print(counttable)
      })
      estims = rbind(estims, data.frame(Riskratio = rr[2,5],
                                        Lower = rr[2,6],
                                        Upper = rr[2,7],
                                        Pval = rr[2, 8],
                                        Feature = feature))
    } else {
      cat('MAF:', maf, 'vartype:', vartype, 'Method:', method, 'Feature:', feature, 'Feature var:', featurevar, '\n')
    }
  }
  ## Include the method and variant type in return data frame for ease of plotting
  estims = estims %>% mutate(Method = method, Type = vartype, Maf = maf)
  return(estims)
}

## Produce all enrichments statistics for a given set of outliers (across all mafs and variant types)
## Input: outlier/control data frame, path to features, whether they are count features, whaether to scale features
## Output: Single data frame with all the enrichments from logistic regression
get.all.enrich <- function(outliers, feature_dir, counts = FALSE, scale = TRUE, suffix = NULL) {
    if (is.null(suffix)) {
        if (counts) {
            feature_dir = paste0(feature_dir, '/counts/')
            suffix = '_counts_bygene.txt'
            scale = FALSE
        } else {
            suffix = '_features_bygene.txt'
        }
    }

    ## get mafs and vartypes - assumes the same vartypes for each MAF
    mafs = list.dirs(feature_dir, full.names = FALSE, recursive = FALSE)
    mafs = c('MAF0-1','MAF1-5','MAF5-10','MAF10-25')
    if (!counts) {
#        mafs = 'MAF0-1'
        mafs = mafs[!(mafs == 'counts')] # removes counts dir, if there
    }
    vartypes = unique(list.dirs(paste0(feature_dir, '/', mafs), full.names = FALSE, recursive = FALSE))
    #vartypes = vartypes[which(vartypes != 'HallLabSV')] # removing for v8
    if (!counts) {
        vartypes = vartypes[which(vartypes != 'HallLabSV')] # Remove SVs for features
    }
    vartypes = c('SNPs','HallLabSV') # adding to combine enrichments for SNPs+indels
    #vartypes = c('SNPs')
    ## for all combinations of mafs and vartypes, merge features with outliers
    ## then collapse into a single data frame
    if (!counts) {
      logit_all = do.call(rbind,
                          mapply(compute.enrichments,
                                 rep(mafs, length(vartypes)), rep(vartypes, each = length(mafs)),
                                 MoreArgs = list(outliers, feature_dir, suffix, scale, counts), SIMPLIFY = FALSE))
    } else {
      logit_all = do.call(rbind,
                          mapply(compute.relative.risk,
                                 rep(mafs, length(vartypes)), rep(vartypes, each = length(mafs)),
                                 MoreArgs = list(outliers, feature_dir, suffix, scale, counts), SIMPLIFY = FALSE))
    }
    return(logit_all)
}

## Function to generate default ggplot colors
gg_color_hue <- function(n) {
    hues = seq(15, 375, length=n+1)
    hcl(h=hues, l=65, c=100)[1:n]
}

## Function for scientific notation 
fancy_scientific <- function(l) {
    ## turn in to character string in scientific notation
    lnew <- format(l, scientific = TRUE)
    ## quote the part before the exponent to keep all the digits
    lnew <- gsub("^(.*)e", "'\\1'e", lnew)
    ## turn the 'e+' into plotmath format
    lnew <- gsub("e", "%*%10^", lnew)
    ## don't use scientific notation if exponent is zero
    for (i in 1:length(lnew)) {
        if (length(grep('+00', lnew[i], fixed = TRUE)) > 0) {
            lnew[i] = l[i]
        }
    }
    ## return this as an expression
    return(parse(text=lnew))
}
