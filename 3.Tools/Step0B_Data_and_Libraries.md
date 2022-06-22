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
<!-- Data -->

```r
onTarget=readRDS('/Users/sinhas8/Project_TrueTarget/Data/onTarget_v2.0.RDS')
KnownTarget_predictions=readRDS('/Users/sinhas8/Project_TrueTarget/Data/KnownTarget_predictions.RDS')
drugCandidates_for_secTargets=readRDS('/Users/sinhas8/Project_TrueTarget/Data/drugCandidates_for_secTargets.RDS')
```
<!-- Functions -->

```r
source('/Users/sinhas8/Project_TrueTarget/Tools/myCustom_functions.R')
```

```
## Warning in file(filename, "r", encoding = encoding): cannot open file '/Users/sinhas8/
## Project_TrueTarget/Tools/myCustom_functions.R': No such file or directory
```

```
## Error in file(filename, "r", encoding = encoding): cannot open the connection
```


Add a new chunk by clicking the *Insert Chunk* button on the toolbar or by pressing *Cmd+Option+I*.

When you save the notebook, an HTML file containing the code and output will be saved alongside it (click the *Preview* button or press *Cmd+Shift+K* to preview the HTML file). 

The preview shows you a rendered HTML copy of the contents of the editor. Consequently, unlike *Knit*, *Preview* does not run any R code chunks. Instead, the output of the chunk when it was last run in the editor is displayed.

