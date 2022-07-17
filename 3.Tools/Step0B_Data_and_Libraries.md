---
title: "R Notebook"
output: html_notebook
---

This is an [R Markdown](http://rmarkdown.rstudio.com) Notebook. When you execute code within the notebook, the results appear beneath the code. 

Try executing this chunk by clicking the *Run* button within the chunk or by placing your cursor inside it and pressing *Cmd+Shift+Enter*. 
<!-- Atomic inputs needed to run a code -->
<!-- Libraries -->

```r
require(parallel); require('tictoc'); library(fgsea); library(ggplot2)
```

```
## Loading required package: parallel
```

```
## Loading required package: tictoc
```

```
## Registered S3 method overwritten by 'data.table':
##   method           from
##   print.data.table
```

```
## Want to understand how all the pieces fit together? Read R for Data Science: https://r4ds.had.co.nz/
```

```r
require(ggpubr); require(statar); require('interactions')
```

```
## Loading required package: ggpubr
```

```
## Loading required package: statar
```

```
## Loading required package: interactions
```

```r
require(stats); require(reshape2); require(pROC)
```

```
## Loading required package: reshape2
```

```r
require(data.table); 
```

```
## Loading required package: data.table
```

```
## data.table 1.14.2 using 1 threads (see ?getDTthreads).  Latest news: r-datatable.com
```

```
## **********
## This installation of data.table has not detected OpenMP support. It should still work but in single-threaded mode.
## This is a Mac. Please read https://mac.r-project.org/openmp/. Please engage with Apple and ask them for support. Check r-datatable.com for updates, and our Mac instructions here: https://github.com/Rdatatable/data.table/wiki/Installation. After several years of many reports of installation problems on Mac, it's time to gingerly point out that there have been no similar problems on Windows or Linux.
## **********
```

```
## 
## Attaching package: 'data.table'
```

```
## The following objects are masked from 'package:reshape2':
## 
##     dcast, melt
```

```r
require(ggrepel)
```

```
## Loading required package: ggrepel
```

```r
require(grid)
```

```
## Loading required package: grid
```

```r
require(cowplot)
```

```
## Loading required package: cowplot
```

```
## 
## Attaching package: 'cowplot'
```

```
## The following object is masked from 'package:ggpubr':
## 
##     get_legend
```

```r
require(gridExtra)
```

```
## Loading required package: gridExtra
```

```r
require('stringb')
```

```
## Loading required package: stringb
```

```r
require(DEGreport)
```

```
## Loading required package: DEGreport
```
<!-- Data -->

```r
setwd('/Users/sinhas8/Project_TrueTarget/Tools_github/')
onTarget=readRDS('../Data/onTarget_v2.0.RDS')
KnownTarget_predictions=readRDS('/Users/sinhas8/Project_TrueTarget/Data/KnownTarget_predictions_v3.RDS')
drugCandidates_for_secTargets=readRDS('../Data/drugCandidates_for_secTargets.RDS')
drugVScrispr_corr_features_list=readRDS('../Data/drugVScrispr_corr_features_list.RDS')
drugVScrispr_corr_Strength=sapply(drugVScrispr_corr_features_list, function(x) x[,2])
corrMat_P=sapply(drugVScrispr_corr_features_list, function(x) x[,1])
corrMat=drugVScrispr_corr_Strength
corrMat_rank=apply(-corrMat, 2, rank)

#Functions needed
source('/Users/sinhas8/Project_TrueTarget/Tools_github/3.Tools/myCustom_functions.R')
```
<!-- Functions for producing negative controld during primary target identification -->

```r
# Create a random df unique to the input df
randomize_DF_with_unique_Pairs<-function(df2randomize=Positive_Set_DrugGene_Pairs){
  df2randomize_collapsed=paste(df2randomize[,1], df2randomize[,2], sep = '_')
  fixedCol=df2randomize[,2]
  shuffledCol=sample(df2randomize[,1])
  new_df=paste(shuffledCol, fixedCol, sep = '_')
  while(sum(df2randomize_collapsed==new_df)>0){
    new_df[df2randomize_collapsed==new_df]=
      paste(sample(df2randomize[,1], sum(df2randomize_collapsed==new_df)), fixedCol[df2randomize_collapsed==new_df], sep = '_')
  }
  df2return=do.call(rbind, strsplit(new_df, '_'))
  df2return=data.frame(df2return[,1], df2return[,2])
  colnames(df2return)=colnames(df2randomize)
  df2return
}

randomize_DF_with_unique_Pairs_V2<-function(df2randomize=Positive_Set_DrugGene_Pairs){
  df2randomize_collapsed=paste(df2randomize[,1], df2randomize[,2], sep = '_')
  fixedCol=df2randomize[,2]
  shuffledCol=sample(rownames(corrMat), nrow(df2randomize))
  new_df=paste(shuffledCol, fixedCol, sep = '_')
  while(sum(df2randomize_collapsed==new_df)>0){
    new_df[df2randomize_collapsed==new_df]=
      paste(sample(rownames(onTarget$corrMat), sum(df2randomize_collapsed==new_df)),
            fixedCol[df2randomize_collapsed==new_df], sep = '_')
  }
  df2return=do.call(rbind, strsplit(new_df, '_'))
  df2return=data.frame(df2return[,1], df2return[,2])
  colnames(df2return)=colnames(df2randomize)
  df2return
}
```
<!-- Compute if a drug is inhibitor or not -->

```r
category_inhibitor=c('inhibitor', 'antagonist', 'blocker')
category_activator=c('agonist', 'activator', 'stimulant')
category_inhibitor_drugs_ID=sapply(category_inhibitor, function(x) grep(x, onTarget$drugCategory$moa))
category_activator_drugs_ID=sapply(category_activator, function(x) grep(x, onTarget$drugCategory$moa))
category_inhibitor_drugs_ID=unique(unlist(category_inhibitor_drugs_ID))
category_activator_drugs_ID=unique(unlist(category_activator_drugs_ID))
onTarget$drugCategory$is.inhibitor=F
onTarget$drugCategory$is.inhibitor[category_inhibitor_drugs_ID]=T
onTarget$drugCategory$is.activator=F
onTarget$drugCategory$is.activator[category_activator_drugs_ID]=T
```

<!-- MATCHED MATRIX -->

```r
matched_cellLines=Reduce(intersect,
                         list(colnames(onTarget$expression_20Q4),
                              colnames(onTarget$avana_22Q2),
                              colnames(onTarget$secondary_prism)))
expression_matched=onTarget$expression_20Q4[,matched_cellLines]
avana_matched=onTarget$avana_22Q2[,matched_cellLines]
drug_matched=onTarget$secondary_prism[,matched_cellLines]
```

Add a new chunk by clicking the *Insert Chunk* button on the toolbar or by pressing *Cmd+Option+I*.

When you save the notebook, an HTML file containing the code and output will be saved alongside it (click the *Preview* button or press *Cmd+Shift+K* to preview the HTML file). 

The preview shows you a rendered HTML copy of the contents of the editor. Consequently, unlike *Knit*, *Preview* does not run any R code chunks. Instead, the output of the chunk when it was last run in the editor is displayed.

