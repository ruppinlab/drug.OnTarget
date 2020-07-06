# Process secondary screen from raw data
require(parallel); require(tictoc); require(WGCNA); require(ggplot2); require(cowplot)
library(tidyr)
setwd('/Users/sinhas8/Project_OffTarget/')
source('drug.OnTarget/3.Tools/myCustom_functions.R')
###########################################################################
# Step 0: Load dataset needed
###########################################################################
secondary=read.csv('secondary-screen-dose-response-curve-parameters.csv')
secondary$col_names=paste(secondary$broad_id, secondary$screen_id, secondary$name,sep='_')

secondary_trimmed=secondary[,c(21, 20, 9)]
secondary_trimmed[c(387225, 387401, 387430, 387436, 387531),]

secondary_trimmed_wide <- spread(secondary_trimmed, col_names, auc)
rownames(secondary_trimmed_wide)=as.character(secondary_trimmed_wide$row_name)
secondary_trimmed_wide=secondary_trimmed_wide[,-1]

column_annotation=t( sapply(colnames(secondary_trimmed_wide), function(x) c(strsplit(x, '_')[[1]])))

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
