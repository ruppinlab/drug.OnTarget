---
title: "Data, libraries and preprocessing"
output: html_notebook
---

<!-- Atomic inputs needed to run a code -->
<!-- Libraries -->

```r
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# 
# BiocManager::install("DEGreport")

### network problem : talked to narisu and talked about the authentication

require(parallel)
require('tictoc')
# library(fgsea)
library(ggplot2)
require(ggpubr); require(statar); require('interactions')
require(stats); require(reshape2); require(pROC)
require(data.table); 
require(ggrepel)
require(grid)
require(cowplot)
require(gridExtra)
require('stringb')
require(DEGreport)
```

```
## Loading required package: DEGreport
```

```
## Warning in library(package, lib.loc = lib.loc, character.only = TRUE, logical.return = TRUE, : there is no package called 'DEGreport'
```

```r
#Functions needed
source('../3.Tools/myCustom_functions.R')
```
<!-- OnTarget Data and the preprocessing-->

```r
# setwd('/Users/sinhanee/Documents/projects/drug.OnTarget/')
onTarget=readRDS('../Data/onTarget_v2.0.RDS')
# KnownTarget_predictions=readRDS('../Data/KnownTarget_predictions_v3.RDS') ## required dataset
# drugCandidates_for_secTargets=readRDS('../Data/drugCandidates_for_secTargets.RDS') 
# drugVScrispr_corr_features_list=readRDS('../Data/drugVScrispr_corr_features_list.RDS')


# drugSubset_map2_annotation <- match(rownames(onTarget$secondary_prism),
#         onTarget$drugCategory$broad_id_trimmed)
# all_drugNames_common = onTarget$drugCategory$name[drugSubset_map2_annotation]
# all_drugTargets = onTarget$drugCategory$target[drugSubset_map2_annotation]
# all_drugMOA = onTarget$drugCategory$moa[drugSubset_map2_annotation]
```
<!-- Run correlation btw all crispr vs all drugs -->

```r
## finding the commom cellines between the CRISPR and drug screen
common_cellLines=intersect(colnames(onTarget$avana_22Q2),colnames(onTarget$secondary_prism))

## extracting names of all drugs and genes
all_genenames=rownames(onTarget$avana_22Q2)
all_drugNames=rownames(onTarget$secondary_prism)

## Running the function to find correlation between all crispr vs all drugs
tic()
drugVScrispr_corr_features_list=mclapply(all_drugNames, 
                                         function(x) correlation_bet_crispr_drug(x), mc.cores = detectCores())
names(drugVScrispr_corr_features_list)=all_drugNames
toc()
```

```
## 441.47 sec elapsed
```

```r
## Save RDS file
# saveRDS(drugVScrispr_corr_features_list,
#         '../Data/drugVScrispr_corr_features_list.RDS')
```

<!-- Get the correlation strength, P value and the rank of the correlation matrix-->

```r
dim(drugVScrispr_corr_Strength)
```

```
## [1] 17386  1618
```

```r
drugVScrispr_corr_Strength=sapply(drugVScrispr_corr_features_list, function(x) x[,2])
corrMat_P=sapply(drugVScrispr_corr_features_list, function(x) x[,1])
corrMat=drugVScrispr_corr_Strength
corrMat_rank=apply(-corrMat, 2, rank)
```
<!-- Preprocessing of Drug Category Data set to identify drug is inhibitor or not -->
<!-- Compute if a drug is inhibitor or not -->

```r
# category_inhibitor=c('inhibitor', 'antagonist', 'blocker')
# category_activator=c('agonist', 'activator', 'stimulant')
# category_inhibitor_drugs_ID=sapply(category_inhibitor, function(x) grep(x, onTarget$drugCategory$moa))
# category_activator_drugs_ID=sapply(category_activator, function(x) grep(x, onTarget$drugCategory$moa))
# category_inhibitor_drugs_ID=unique(unlist(category_inhibitor_drugs_ID))
# category_activator_drugs_ID=unique(unlist(category_activator_drugs_ID))
# onTarget$drugCategory$is.inhibitor=F
# onTarget$drugCategory$is.inhibitor[category_inhibitor_drugs_ID]=T
# onTarget$drugCategory$is.activator=F
# onTarget$drugCategory$is.activator[category_activator_drugs_ID]=T
match(grep('T', onTarget$drugCategory$is.activator), grep('T', drug_cat$is.activator.check)) ## tell the activator part to sanju is.activator ~1640, is.activator.check ~800.
```

```
## integer(0)
```

```r
drug_cat = onTarget$drugCategory %>% mutate(is.inhibitor.check = F) 
drug_cat= drug_cat %>% dplyr::mutate(is.inhibitor.check = ifelse(grepl('inhibitor', drug_cat$moa) |
                                                                   grepl('antagonist', drug_cat$moa) |
                                                                   grepl('blocker', drug_cat$moa) ,TRUE, FALSE),
                                     is.activator.check = ifelse(grepl('activator', drug_cat$moa) |
                                                                   grepl('\\bagonist\\b', drug_cat$moa) |
                                                                   grepl('stimulant', drug_cat$moa) ,TRUE, FALSE))
```

<!-- MATCHED MATRIX for expression data -->

```r
matched_cellLines=Reduce(intersect,
                         list(colnames(onTarget$expression_20Q4),
                              colnames(onTarget$avana_22Q2),
                              colnames(onTarget$secondary_prism)))
expression_matched=onTarget$expression_20Q4[,matched_cellLines]
avana_matched=onTarget$avana_22Q2[,matched_cellLines]
drug_matched=onTarget$secondary_prism[,matched_cellLines]
```

