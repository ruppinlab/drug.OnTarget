# Landscape of shRNA
###########################################################################
# Step 1B: Target Correlation
###########################################################################
Targets_list=lapply(as.character(onTarget$Annotated_Target),
                    function(x) gsub(' ','',strsplit(x, '\\,')[[1]]) )
Annotated_Target_corr=sapply(1:length(Targets_list), function(x)
  sapply(Targets_list[[x]], function(y)  err_handle(onTarget$corrMat_shRNA[y,x]) ))
Annotated_Target_corr_mean=sapply(Annotated_Target_corr, function(x) mean(x, na.rm=T))
Annotated_Target_corr_max=sapply(Annotated_Target_corr, function(x) max(x, na.rm=T))
tiff('/Users/sinhas8/Project_OffTarget/4.Results/Target_annotation_corMat_landscape_shRNA.tiff', 
     height=600, width = 600)
plot_grid(plotHist(Annotated_Target_corr_mean, xtitle='Annotated_Target_corr_mean'),
          plotHist(Annotated_Target_corr_max,xtitle='Annotated_Target_corr_max'),
          plotHist(Annotated_Target_corr_rank_mean, xtitle='Annotated_Target_corr_mean_RANK'),
          plotHist(Annotated_Target_corr_rank_max,xtitle='Annotated_Target_corr_max_RANK'))
dev.off()
###########################################################################
# Step 1C: Target Corr Rank 
###########################################################################
Annotated_Target_corr_rank=sapply(1:length(Targets_list), function(x)
  sapply(Targets_list[[x]], function(y)  err_handle(onTarget$corrMat_shRNA_rank[y,x]) ))
Annotated_Target_corr_rank_mean=sapply(Annotated_Target_corr_rank, function(x) mean(x, na.rm=T))
Annotated_Target_corr_rank_max=sapply(Annotated_Target_corr_rank, function(x) max(x, na.rm=T))
###########################################################################
# Perfect annotation
###########################################################################
Annotated_Target_corr[which(Annotated_Target_corr_rank_mean== max(Annotated_Target_corr_rank_mean, na.rm = T) )]
# Semi-Perfect annotation
Annotated_Target_corr[which(Annotated_Target_corr_rank_mean>max(Annotated_Target_corr_rank_mean, na.rm = T)-1000)]
# Semi-Perfect annotation
(Annotated_Target_corr[depmap$drug_target_map$phase=='Launched'])[which(Annotated_Target_corr_rank_mean[depmap$drug_target_map$phase=='Launched']>max(Annotated_Target_corr_rank_mean, na.rm = T)-200)]
###########################################################################
# Recent Remapping
###########################################################################
# SGI-1776
cmat_rank['PIM1',grep('1776', depmap$drug_target_map$name, ignore.case = T)]
# nutlin-3a to TP53
which.max(cmat_rank[,grep('nutlin-3', depmap$drug_target_map$name, ignore.case = T)])

