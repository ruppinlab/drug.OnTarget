# Process secondary screen from raw data
require(parallel); require(tictoc); require(WGCNA); require(ggplot2); require(cowplot); library(tidyr)
setwd('/Users/sinhas8/Project_OffTarget/')
source('drug.OnTarget/3.Tools/myCustom_functions.R')
###########################################################################
# Step 0: Load dataset needed
###########################################################################
secondary=read.csv('/Users/sinhas8/Project_OffTarget/secondary-screen-dose-response-curve-parameters.csv')
secondary$col_names=paste(secondary$broad_id, secondary$screen_id, secondary$name, secondary$phase, sep='_')
head(secondary)
# Getting AUC
secondary_trimmed=secondary[,c(21, 20, 9)]
secondary_trimmed_wide <- spread(secondary_trimmed, col_names, auc)
rownames(secondary_trimmed_wide)=as.character(secondary_trimmed_wide$row_name)
secondary_trimmed_wide=secondary_trimmed_wide[,-1]

column_annotation=t(sapply(colnames(secondary_trimmed_wide), function(x) c(strsplit(x, '_')[[1]])))

column_annotation=data.frame(column_annotation)
colnames(column_annotation)=c('Broad_id', 'Screen_id', 'CommonName')
rownames(column_annotation)=NULL
column_annotation$Broad_id_trimmed=sapply(as.character(column_annotation$Broad_id), function(x) strsplit(x, '-')[[1]][2])
secondary_trimmed_wide_t=t(secondary_trimmed_wide)
rownames(secondary_trimmed_wide_t)=column_annotation$Broad_id_trimmed
secondary=list(secondary_screen=secondary_trimmed_wide_t,
               drug_Info=column_annotation,
               cellLine_Info=rownames(secondary_trimmed_wide))
saveRDS(secondary, '2.Data/secondary_screen_processed.RDS')
###########################################################################
# Extract IC50
###########################################################################
secondary=read.csv('/Users/sinhas8/Project_OffTarget/secondary-screen-dose-response-curve-parameters.csv')
secondary$col_names=paste(secondary$broad_id, secondary$screen_id, secondary$name, secondary$phase, sep='_')
head(secondary)
secondary_trimmed=secondary[,c(21, 20, 11)]
secondary_trimmed_wide <- spread(secondary_trimmed, col_names, ic50)
rownames(secondary_trimmed_wide)=as.character(secondary_trimmed_wide$row_name)
secondary_trimmed_wide=secondary_trimmed_wide[,-1]

column_annotation=t(sapply(colnames(secondary_trimmed_wide), function(x) c(strsplit(x, '_')[[1]])))

column_annotation=data.frame(column_annotation)
colnames(column_annotation)=c('Broad_id', 'Screen_id', 'CommonName')
rownames(column_annotation)=NULL
column_annotation$Broad_id_trimmed=sapply(as.character(column_annotation$Broad_id), function(x) strsplit(x, '-')[[1]][2])
secondary_trimmed_wide_t=t(secondary_trimmed_wide)
rownames(secondary_trimmed_wide_t)=column_annotation$Broad_id_trimmed
onTarget$secondary_prism_ic50=secondary_trimmed_wide_t

