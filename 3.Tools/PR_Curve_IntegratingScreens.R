# Precision-Recall for integrated genetic screen based predictions
# Maximizing a Precision recall curve given a gold-standard (Hits with citation>3)
###########################################################################
# Add Cripsr & shRNA corr Score for known Targets
###########################################################################
require(ggpubr)
require(pROC)
source('drug.OnTarget/3.Tools/drug_type_effect.R')
df2plot$both_corr=sapply(1:nrow(df2plot), function(x)
  err_handle(onTarget$corrMat_bothScreens[as.character(df2plot$Drugbank_Gene)[x], as.character(df2plot$drugName)[x]]) ) 
df2plot$both_corr_rank=sapply(1:nrow(df2plot), function(x) 
  err_handle(onTarget$corrMat_rank[as.character(df2plot$Drugbank_Gene)[x], as.character(df2plot$drugName)[x]]) ) 
###########################################################################
# Subset showing signal from DrugBank
###########################################################################
test_AUC<-function(Citations_Count_thr=4, seedNumber=1, Count_KnownTarget_thr=3){
  cond1 = df2plot$Action %in% c('inhibitor', 'inhibitor, competitive',
                                'antagonist', 'negative modulator')
  cond2 = df2plot$drugCategory %in% 'targeted cancer'
  cond3=df2plot$Citations_Count>Citations_Count_thr
  cond4=df2plot$Count_KnownTarget<Count_KnownTarget_thr
  sum(cond1 & cond2 & cond3 & cond4)
  drugs_forPR_try1=df2plot[cond1 & cond2 & cond3 & cond4,]
  Positive_Set_DrugGene_Pairs=drugs_forPR_try1[,1:2]
  set.seed(seedNumber)
  Negative_Set_DrugGene_Pairs=randomize_DF_with_unique_Pairs(Positive_Set_DrugGene_Pairs)
  Positive_Set_DrugGene_Pairs$Label=1
  Negative_Set_DrugGene_Pairs$Label=0
  Complete_Set=rbind(Positive_Set_DrugGene_Pairs, Negative_Set_DrugGene_Pairs)
  
  Complete_Set$drug_BroadID=onTarget$drugCategory$broad_id_trimmed[match(tolower(Complete_Set$Chemical.Name), tolower(onTarget$drugCategory$name))]
  Complete_Set$Label= factor(Complete_Set$Label, labels = c('Negative', 'Positive'))
  Complete_Set$both_corr=sapply(1:nrow(Complete_Set), function(x) 
    err_handle(onTarget$corrMat_bothScreens_rank[as.character(Complete_Set$Drugbank_Gene)[x], as.character(Complete_Set$drug_BroadID)[x]]))
  roc_curve=roc(Complete_Set$Label, Complete_Set$both_corr)$auc
  roc_curve
}
sapply(2:5, function(x) test_AUC(Citations_Count_thr=x, Count_KnownTarget_thr=2))
# tiff('/Users/sinhas8/Project_OffTarget/4.Results/crisprVSshRNA_prediction_Figure3A.tiff')
# ggplot(df2plot_forPlot, aes(y=Corr_Rho, fill=factor(Label), x=factor(Label)))+
#   geom_boxplot()+
#   facet_wrap(~Screen)+
#   stat_compare_means(method='wilcox')+
#   theme_bw(base_size = 25)+
#   labs(y='Correlation Rho', x='Ground Truth', fill='Ground Truth')+
#   theme(legend.position = 'top')
# dev.off()
# 
# A<-ggroc(roc(crispr_score$Label, crispr_score$Corr_Rho))+
#   theme_bw(base_size = 20)+
#   ggtitle('CRISPR-Based; AUC=0.73')
# B<-ggroc(roc(shRNA_score$Label, shRNA_score$Corr_Rho))+
#   theme_bw(base_size = 20)+
#   ggtitle('shRNA-Based; AUC=0.68')
# tiff('/Users/sinhas8/Project_OffTarget/4.Results/crisprVSshRNA_AUC_Figure3B.tiff')
# grid.arrange(A, B)
# dev.off()
