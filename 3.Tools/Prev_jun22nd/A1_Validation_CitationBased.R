# A1 gold-standard (Hits with citation>3)
# ValidataionA1 - By Curating BioGrid Chemical Interaction database we would like to test 
require(stats)
require(reshape2)
require(statar)
require(ggplot2)
require(ggpubr)
require(pROC)
biogrid=read.csv('/Users/sinhas8/Project_OffTarget/2.Data/BIOGRID-CHEMICALS-3.5.182.chemtab.txt',
                 sep='\t')
common_Drug_lowerCase=tolower(levels(biogrid$Chemical.Name))[
  !is.na(match(tolower(levels(biogrid$Chemical.Name)), tolower(onTarget$PredvsKnown_scores$CommonDrugName) ))]
biogrid_matched=biogrid[!is.na(match(tolower(biogrid$Chemical.Name),
                                     common_Drug_lowerCase)),]
PredvsKnown_scores_matched=onTarget$PredvsKnown_scores[
  !is.na(match(tolower(onTarget$PredvsKnown_scores$CommonDrugName), common_Drug_lowerCase)),]
biogrid_matched$Chemical.Name= factor(as.character(biogrid_matched$Chemical.Name))
biogrid_matched_Subsetlist=split(biogrid_matched[,c('Official.Symbol', 'Chemical.Name',
                                                    'Interaction.Type', 'Action')],
                                 biogrid_matched$Chemical.Name)
biogrid_matched_Subsetlist_reformat=lapply(biogrid_matched_Subsetlist, function(x)
  dcast(x, Official.Symbol + Chemical.Name + Action ~  Interaction.Type))

# dcast(biogrid_matched_Subsetlist$Barbital, Official.Symbol + Chemical.Name + Action ~  Interaction.Type)
biogrid_matched_Subsetlist_withScore=do.call(rbind, biogrid_matched_Subsetlist_reformat)
colnames(biogrid_matched_Subsetlist_withScore)[4]='Score'
biogrid_matched_Subsetlist_withScore$Score=as.numeric(biogrid_matched_Subsetlist_withScore$Score)
biogrid_matched_Subsetlist_withScore$Score[is.na(biogrid_matched_Subsetlist_withScore$Score)]=1
biogrid_matched_Subsetlist_withScore=biogrid_matched_Subsetlist_withScore[
  order(biogrid_matched_Subsetlist_withScore$Score, decreasing = T),]
###########################################################################
# Step 2
###########################################################################
DrugBankScore_vsOurScore=data.frame(biogrid_matched_Subsetlist_withScore, 
                                    PredvsKnown_scores_matched[match(
                                      tolower(biogrid_matched_Subsetlist_withScore$Chemical.Name),
                                      tolower(PredvsKnown_scores_matched$CommonDrugName)),])

DrugBankScore_vsOurScore$drugCategory=onTarget$drugCategory$drug_category[
  match(DrugBankScore_vsOurScore$drugName, colnames(onTarget$corrMat))]

colnames(DrugBankScore_vsOurScore)[c(1)]='Drugbank_Gene'
colnames(DrugBankScore_vsOurScore)[c(4)]='Citations_Count'

df2plot=as.data.frame(DrugBankScore_vsOurScore)
df2plot$Count_KnownTarget=sapply(df2plot$KnownTarget,
                                 function(x) length(strsplit(as.character(x), '\\,')[[1]])) 
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
  print(dim(Positive_Set_DrugGene_Pairs))
  set.seed(seedNumber)
  Negative_Set_DrugGene_Pairs=randomize_DF_with_unique_Pairs(Positive_Set_DrugGene_Pairs)
  Positive_Set_DrugGene_Pairs$Label=1
  Negative_Set_DrugGene_Pairs$Label=0
  Complete_Set=rbind(Positive_Set_DrugGene_Pairs, Negative_Set_DrugGene_Pairs)
  
  Complete_Set$drug_BroadID=onTarget$drugCategory$broad_id_trimmed[match(tolower(Complete_Set$Chemical.Name), tolower(onTarget$drugCategory$name))]
  Complete_Set$Label= factor(Complete_Set$Label, labels = c('Negative', 'Positive'))
  # Complete_Set$both_corr=sapply(1:nrow(Complete_Set), function(x) 
  #   err_handle(onTarget$corrMat_bothScreens_rank[as.character(Complete_Set$Drugbank_Gene)[x], as.character(Complete_Set$drug_BroadID)[x]]))
  # roc_curve=roc(Complete_Set$Label, Complete_Set$both_corr)$auc
  # roc_curve
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
sapply(2:7, function(x) test_AUC(Citations_Count_thr=x, Count_KnownTarget_thr=2))
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
