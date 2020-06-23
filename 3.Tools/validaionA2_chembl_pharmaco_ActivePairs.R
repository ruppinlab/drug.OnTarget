# Validation V2 - Gold Standard based on pharmacologically active drugs from DrugBank
require(ggpubr)
load('/Users/sinhas8/Downloads/drugs2targets.RData')

numberOfTargets=table(drugs2targets$drug.name)
drugs2targets$numberOfTargets=numberOfTargets[match(drugs2targets$drug.name, names(numberOfTargets))]

active_drugs2targets=drugs2targets[which(drugs2targets$pharmaco.active=='yes' & 
                                           drugs2targets$action.simp.strict=='inhibition'),]
numberOfActiveTargets=table(as.character(active_drugs2targets$drug.name))
active_drugs2targets$numberOfActiveTargets=numberOfActiveTargets[match(active_drugs2targets$drug.name, names(numberOfActiveTargets))]
matched_gs3=active_drugs2targets[(tolower(active_drugs2targets$drug.name) %in%  tolower(onTarget$drugCategory$name)),]
matched_gs3=na.omit(matched_gs3)
head(matched_gs3)

# The AUC is super low :p 
test_AUC<-function(seedNumber=1, numberOfActiveTargets_Thr=2, numberOfTargets_Thr=2){
  cond1=matched_gs3$numberOfActiveTargets<numberOfActiveTargets_Thr
  cond2=matched_gs3$numberOfTargets<numberOfTargets_Thr
  
  Positive_Set_DrugGene_Pairs=na.omit(data.frame(Drugbank_Gene=matched_gs3$protein.gene.symbol[cond1 & cond2], 
                                                 Chemical.Name=matched_gs3$drug.name[cond1 & cond2]))
  Positive_Set_DrugGene_Pairs=unique(Positive_Set_DrugGene_Pairs)
  Positive_Set_DrugGene_Pairs=Positive_Set_DrugGene_Pairs[Positive_Set_DrugGene_Pairs$Drugbank_Gene %in%
                                rownames(onTarget$corrMat_bothScreens),]
  print(dim(Positive_Set_DrugGene_Pairs))
  set.seed(seedNumber)
  Negative_Set_DrugGene_Pairs=randomize_DF_with_unique_Pairs(Positive_Set_DrugGene_Pairs)
  Positive_Set_DrugGene_Pairs$Label=1
  Negative_Set_DrugGene_Pairs$Label=0
  Complete_Set=rbind(Positive_Set_DrugGene_Pairs, Negative_Set_DrugGene_Pairs)
  
  Complete_Set$Label= factor(Complete_Set$Label, labels = c('Negative', 'Positive'))
  Complete_Set$drug_BroadID=onTarget$drugCategory$broad_id_trimmed[match(tolower(Complete_Set$Chemical.Name), tolower(onTarget$drugCategory$name))]
  Complete_Set$crispr_corr=sapply(1:nrow(Complete_Set), function(x)
    err_handle(onTarget$corrMat[as.character(Complete_Set$Drugbank_Gene)[x], as.character(Complete_Set$drug_BroadID)[x]]))
  Complete_Set$shRNA_corr=sapply(1:nrow(Complete_Set), function(x)
    err_handle(onTarget$corrMat_shRNA[as.character(Complete_Set$Drugbank_Gene)[x], as.character(Complete_Set$drug_BroadID)[x]]))
  Complete_Set$both_corr=sapply(1:nrow(Complete_Set), function(x)
    err_handle(onTarget$corrMat_bothScreens[as.character(Complete_Set$Drugbank_Gene)[x], as.character(Complete_Set$drug_BroadID)[x]]))
  Complete_Set=na.omit(Complete_Set)
  roc_curve_crispr=roc(Complete_Set$Label, Complete_Set$crispr_corr)$auc
  roc_curve_shRNA=roc(Complete_Set$Label, Complete_Set$shRNA_corr)$auc
  roc_curve_both=roc(Complete_Set$Label, Complete_Set$both_corr)$auc
  c(roc_curve_crispr, roc_curve_shRNA, roc_curve_both)
}
test_AUC(seedNumber=2, numberOfTargets_Thr=2, numberOfActiveTargets_Thr=2)

###########################################################################
# Figure 3A - Gold-set scores
###########################################################################
# crispr_score=data.frame(Complete_Set[,-6], Screen='CRISPR')
# colnames(crispr_score)[5]='Corr_Rho'
# shRNA_score=data.frame(Complete_Set[,-5], Screen='shRNA')
# colnames(shRNA_score)[5]='Corr_Rho'
# df2plot_forPlot=rbind(crispr_score, shRNA_score)
# df2plot_forPlot$Label= factor(df2plot_forPlot$Label, labels = c('Negative', 'Positive'))

# tiff('/Users/sinhas8/Project_OffTarget/4.Results/crisprVSshRNA_prediction_Figure3A_PharmaActive_DrugBank.tiff')
# ggplot(df2plot_forPlot, aes(y=Corr_Rho, fill=factor(Label), x=factor(Label)))+
#   geom_boxplot()+
#   facet_wrap(~Screen)+
#   stat_compare_means(method='wilcox')+
#   theme_bw(base_size = 25)+
#   labs(y='Correlation Rho', x='Ground Truth', fill='Ground Truth')+
#   theme(legend.position = 'top')
# dev.off()
# 
# require(pROC)
# A<-ggroc(roc(crispr_score$Label, crispr_score$Corr_Rho))+
#   theme_bw(base_size = 20)+
#   ggtitle('CRISPR-Based; AUC=0.66')
# B<-ggroc(roc(shRNA_score$Label, shRNA_score$Corr_Rho))+
#   theme_bw(base_size = 20)+
#   ggtitle('shRNA-Based; AUC=0.67')
# tiff('/Users/sinhas8/Project_OffTarget/4.Results/crisprVSshRNA_AUC_Figure3B_PharmaActiveCompounds.tiff')
# grid.arrange(A, B)
# dev.off()
