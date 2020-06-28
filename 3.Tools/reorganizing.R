source('/Users/sinhas8/myCustom_functions.R')
##################################################
# Reorganizing drug cor data
##################################################
onTarget=readRDS('/Users/sinhas8/Project_OffTarget/2.Data/onTarget.RDS')
onTarget=list()
onTarget$corrMat=cmat
onTarget$drugsCommonName=depmap$drug_target_map$name
onTarget$Annotated_Target=depmap$drug_target_map$target
onTarget$broad_id=depmap$drug_target_map$broad_id
onTarget$expression=depmap$expression
onTarget$mutations_matrix=depmap$mutations_matrix
onTarget$avana=depmap$avana
onTarget$achilles=depmap$achilles
onTarget$drug_prism=depmap$drug_prism
onTarget$annotation=depmap$annotation
onTarget$drugCategory=read.csv('/Users/sinhas8/Project_OffTarget/2.Data/drugCategory.csv')
onTarget$drugCategory=onTarget$drugCategory[match(colnames(onTarget$corrMat), 
                            sapply(as.character(onTarget$drugCategory$broad_id), function(x) strsplit(x, '-')[[1]][2])),]
##################################################
# Predicted best Target and Drug - CRISPR
##################################################
onTarget$corrMat_rank=apply(onTarget$corrMat, 2, rank)
predicted_Target = apply(onTarget$corrMat, 2,
                         function(x) c(PredTarget=rownames(onTarget$corrMat)[which(x==max(x))], 
                                       Score=max(x)))
onTarget$Top_predicted_Target=data.frame(drugName=colnames(onTarget$corrMat),
                                         t(predicted_Target))
Top_predicted_Drug = apply(onTarget$corrMat, 1,
                         function(x) c(PredTarget=colnames(onTarget$corrMat)[which(x==max(x))], 
                                       Score=max(x)))
onTarget$Top_predicted_Drug = data.frame(GeneName=rownames(onTarget$corrMat),
                                         t(Top_predicted_Drug))
# saveRDS(onTarget, '/Users/sinhas8/Project_OffTarget/2.Data/onTarget.RDS')
##################################################
#Phase 2: Add for shRNa:
##################################################
cmat=readRDS('/Users/sinhas8/Project_OffTarget/2.Data/corrMat_shRNA.RDS')
onTarget$corrMat_shRNA = cmat
colnames(onTarget$corrMat_shRNA)=colnames(onTarget$corrMat)
rownames(onTarget$corrMat_shRNA)=rownames(onTarget$achilles)
onTarget$corrMat_shRNA_rank=apply(onTarget$corrMat_shRNA, 2, function(x) rank(x) )
onTarget$drugCategory$broad_id_trimmed=colnames(onTarget$corrMat)
##################################################
# Predicted best Target and Drug - shRNA
##################################################
predicted_Target_shRNA = apply(onTarget$corrMat_shRNA, 2,
                         function(x) c(PredTarget=rownames(onTarget$corrMat_shRNA)[which(x==max(x))], 
                                       Score=max(x)))
onTarget$Top_predicted_Target_shRNA=data.frame(drugName=colnames(onTarget$corrMat_shRNA),
                                         t(predicted_Target_shRNA))
Top_predicted_Drug_shRNA = apply(onTarget$corrMat_shRNA, 1,
                           function(x) c(PredTarget=colnames(onTarget$corrMat_shRNA)[which(x==max(x))], 
                                         Score=max(x)))
onTarget$Top_predicted_Drug_shRNA = data.frame(GeneName=rownames(onTarget$corrMat_shRNA),
                                         t(Top_predicted_Drug_shRNA))


# Add another variable onTarget
onTarget$PredvsKnown_scores=cbind(CommonDrugName=onTarget$drugsCommonName[match(df2plot$drugName, colnames(onTarget$corrMat))], df2plot)
colnames(onTarget$PredvsKnown_scores)[5:10]=c('KnownTarget', 'KnownTarget_corrMean', 'KnownTarget_corrMax',
                                              'KnownTarget_corrRank_mean', 'KnownTarget_corrRank_min',
                                              'Best_among_KnownTarget_based_onCorr')
##################################################
# Annotated_Target_corrMeasures - CRISPR
##################################################
Targets_list = lapply(as.character(onTarget$Annotated_Target),
                      function(x) gsub(' ','',strsplit(x, '\\,')[[1]]) )
Annotated_Target_corr=sapply(1:length(Targets_list), function(x)
  sapply(Targets_list[[x]], function(y)  err_handle(onTarget$corrMat[y,x]) ))
Annotated_Target_corr_mean=sapply(Annotated_Target_corr, function(x) mean(x, na.rm=T))
Annotated_Target_corr_max=sapply(Annotated_Target_corr, function(x) max(x, na.rm=T))
Annotated_Target_corr_Rank=sapply(1:length(Targets_list), function(x)
  sapply(Targets_list[[x]], function(y)  err_handle(onTarget$corrMat_rank[y,x]) ))
Annotated_Target_corrRank_mean=sapply(Annotated_Target_corr_Rank, function(x) mean(x, na.rm=T))
Annotated_Target_corrRank_min=18334-sapply(Annotated_Target_corr_Rank, function(x) max(x, na.rm=T))
Best_among_Annotated_Target=sapply(Annotated_Target_corr, function(x) 
  if(length(which.max(x))==0){NA} else names(x)[which.max(x)] )
onTarget$Annotated_Target_corrMeasures=data.frame(Annotated_Target=onTarget$Annotated_Target,
                                                  corr_mean=Annotated_Target_corr_mean,
                                                  corr_max=Annotated_Target_corr_max,
                                                  corrRank_mean=Annotated_Target_corr_mean,
                                                  corrRank_min=Annotated_Target_corrRank_min,
                                                  Best_among_Annotated_Target,
                                                  drugCategory=onTarget$drugCategory$drug_category)
###########################################################################
# Annotated_Target_corrMeasures - shRNA
###########################################################################
Annotated_Target_corr=sapply(1:length(Targets_list), function(x)
  sapply(Targets_list[[x]], function(y)  err_handle(onTarget$corrMat_shRNA[y,x]) ))
Annotated_Target_corr_mean=sapply(Annotated_Target_corr, function(x) mean(x, na.rm=T))
Annotated_Target_corr_max=sapply(Annotated_Target_corr, function(x) max(x, na.rm=T))
Annotated_Target_corr_Rank=sapply(1:length(Targets_list), function(x)
  sapply(Targets_list[[x]], function(y)  err_handle(onTarget$corrMat_shRNA_rank[y,x]) ))
Annotated_Target_corrRank_mean=sapply(Annotated_Target_corr_Rank, function(x) mean(x, na.rm=T))
Annotated_Target_corrRank_min=nrow(onTarget$corrMat_shRNA)+1 - sapply(Annotated_Target_corr_Rank, function(x) max(x, na.rm=T))
Best_among_Annotated_Target=sapply(Annotated_Target_corr, function(x) 
  if(length(which.max(x))==0){NA} else names(x)[which.max(x)] )
onTarget$Annotated_Target_corrMeasures_shRNA=data.frame(Annotated_Target=onTarget$Annotated_Target,
                                                        corr_mean=Annotated_Target_corr_mean,
                                                        corr_max=Annotated_Target_corr_max,
                                                        corrRank_mean=Annotated_Target_corr_mean,
                                                        corrRank_min=Annotated_Target_corrRank_min,
                                                        Best_among_Annotated_Target,
                                                        drugCategory=onTarget$drugCategory$drug_category)

###########################################################################
# Version2
###########################################################################
saveRDS(onTarget, '/Users/sinhas8/Project_OffTarget/2.Data/onTarget_v2.RDS')
###########################################################################
# Version3-Changes
###########################################################################
# Integrating shRNA and CRISPR genetic screen to create a single score (max from both)
setwd('/Users/sinhas8/Project_OffTarget/')
source('drug.OnTarget/3.Tools/Root_Step0.R')
source('drug.OnTarget/3.Tools/myCustom_functions.R')
dim(onTarget$corrMat); dim(onTarget$corrMat_shRNA)
# Unequal number of genes; For Genes common to both
# we will take the max corr out of two screen (RNAi and CRISPR)
common_genes=intersect(rownames(onTarget$corrMat), rownames(onTarget$corrMat_shRNA))
crispr_unique_genes=rownames(onTarget$corrMat)[!(rownames(onTarget$corrMat) %in% rownames(onTarget$corrMat_shRNA))]
shRNA_unique_genes=rownames(onTarget$corrMat_shRNA)[!(rownames(onTarget$corrMat_shRNA) %in% rownames(onTarget$corrMat))]

infunc_mat=rbind(onTarget$corrMat[crispr_unique_genes,],onTarget$corrMat_shRNA[shRNA_unique_genes,])
corrMat_integrated=sapply(common_genes, function(x) pmax(onTarget$corrMat[x,], onTarget$corrMat_shRNA[x,], na.rm = T))
corrMat_integrated_t=t(corrMat_integrated)
corrMat_bothScreens=rbind(corrMat_integrated_t, infunc_mat)
corrMat_bothScreens=corrMat_bothScreens[order(rownames(corrMat_bothScreens)),]
onTarget$corrMat_bothScreens=corrMat_bothScreens
onTarget$corrMat_shRNA_rank=apply(onTarget$corrMat_shRNA, 2, function(x) rank(x) )
onTarget$drugCategory$broad_id_trimmed=colnames(onTarget$corrMat)
###########################################################################
# Version3-Saving
###########################################################################
saveRDS(onTarget, '/Users/sinhas8/Project_OffTarget/2.Data/onTarget_v3.RDS')

# add the rank matrix
onTarget$corrMat_bothScreens_rank=apply(onTarget$corrMat_bothScreens, 2, function(x) rank(x) )
##################################################
# Predicted best Target and Drug - Integrated Both
##################################################
predicted_Target_bothScreens = apply(onTarget$corrMat_bothScreens, 2,
                               function(x) c(PredTarget=rownames(onTarget$corrMat_bothScreens)[which(x==max(x))], 
                                             Score=max(x)))
max(c(0.1, NA))
onTarget$Top_predicted_Target_bothScreens=data.frame(drugName=colnames(onTarget$corrMat_bothScreens),
                                               t(predicted_Target_bothScreens))
Top_predicted_Drug_bothScreens = apply(onTarget$corrMat_bothScreens, 1,
                                 function(x) c(PredTarget=colnames(onTarget$corrMat_bothScreens)[which(x==max(x))], 
                                               Score=max(x)))
onTarget$Top_predicted_Drug_bothScreens = data.frame(GeneName=rownames(onTarget$corrMat_bothScreens),
                                               t(Top_predicted_Drug_bothScreens))

# Add another variable onTarget
# onTarget$PredvsKnown_scores=cbind(CommonDrugName=onTarget$drugsCommonName[match(df2plot$drugName, colnames(onTarget$corrMat))], df2plot)
# colnames(onTarget$PredvsKnown_scores)[5:10]=c('KnownTarget', 'KnownTarget_corrMean', 'KnownTarget_corrMax',
#                                               'KnownTarget_corrRank_mean', 'KnownTarget_corrRank_min',
#                                               'Best_among_KnownTarget_based_onCorr')

###########################################################################
# Annotated_Target_corrMeasures - bothScreens
###########################################################################
Targets_list=sapply(onTarget$drugCategory$target, function(x) 
  strsplit(as.character(x), ', '))
Annotated_Target_corr=sapply(1:length(Targets_list), function(x)
  sapply(Targets_list[[x]], function(y)  err_handle(onTarget$corrMat_bothScreens[y,x]) ))
Annotated_Target_corr_mean=sapply(Annotated_Target_corr, function(x) mean(x, na.rm=T))
Annotated_Target_corr_max=sapply(Annotated_Target_corr, function(x) max(x, na.rm=T))
Annotated_Target_corr_Rank=sapply(1:length(Targets_list), function(x)
  sapply(Targets_list[[x]], function(y)  err_handle(onTarget$corrMat_bothScreens_rank[y,x]) ))
Annotated_Target_corrRank_mean=sapply(Annotated_Target_corr_Rank, function(x) mean(x, na.rm=T))
Annotated_Target_corrRank_min=nrow(onTarget$corrMat_bothScreens)+1 - sapply(Annotated_Target_corr_Rank, function(x) max(x, na.rm=T))
Best_among_Annotated_Target=sapply(Annotated_Target_corr, function(x) 
  if(length(which.max(x))==0){NA} else names(x)[which.max(x)] )
onTarget$Annotated_Target_corrMeasures_bothScreens=data.frame(Annotated_Target=onTarget$Annotated_Target,
                                                        corr_mean=Annotated_Target_corr_mean,
                                                        corr_max=Annotated_Target_corr_max,
                                                        corrRank_mean=Annotated_Target_corr_mean,
                                                        corrRank_min=Annotated_Target_corrRank_min,
                                                        Best_among_Annotated_Target,
                                                        drugCategory=onTarget$drugCategory$drug_category)

###########################################################################
# Add PredvsKnown for integrated Predicting
###########################################################################
onTarget$PredvsKnown_scores_bothScreens=onTarget$Annotated_Target_corrMeasures_bothScreens
colnames(onTarget$PredvsKnown_scores_bothScreens)=c('KnownTarget', 'KnownTarget_corrMean', 'KnownTarget_corrMax',
                                                    'KnownTarget_corrRank_mean', 'KnownTarget_corrRank_min',
                                                    'Best_among_KnownTarget_based_onCorr', 'drugCategory')
onTarget$PredvsKnown_scores_bothScreens=data.frame(CommonDrugName=onTarget$drugCategory$name,
                                                   onTarget$Top_predicted_Target_bothScreens,
                                                   onTarget$PredvsKnown_scores_bothScreens)
onTarget$PredvsKnown_scores_bothScreens$Score=as.numeric(as.character(onTarget$PredvsKnown_scores_bothScreens$Score))
###########################################################################
# Adding z-score
###########################################################################
onTarget$corrMat_zscore=apply(onTarget$corrMat, 1, scale)
onTarget$corrMat_zscore=t(onTarget$corrMat_zscore)
colnames(onTarget$corrMat_zscore)=colnames(onTarget$corrMat)
onTarget$corrMat_shRNA_zscore=apply(onTarget$corrMat_shRNA, 1, scale)
onTarget$corrMat_shRNA_zscore=t(onTarget$corrMat_shRNA_zscore)
colnames(onTarget$corrMat_shRNA_zscore)=colnames(onTarget$corrMat_shRNA)

onTarget$corrMat_bothScreens_zscore=apply(onTarget$corrMat_bothScreens, 1, scale)
onTarget$corrMat_bothScreens_zscore=t(onTarget$corrMat_bothScreens_zscore)
colnames(onTarget$corrMat_bothScreens_zscore)=colnames(onTarget$corrMat_bothScreens)

###########################################################################
# Version3-Saving
###########################################################################
saveRDS(onTarget, '/Users/sinhas8/Project_OffTarget/2.Data/onTarget_v3.RDS')
###########################################################################
# Version 4- Creating
###########################################################################


