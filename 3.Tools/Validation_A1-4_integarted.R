# Integrating A1-A4
# A4
seedNumber=1; numberOfTargets_Thr=2; Avg_Rating_thr=3
cond1=sab_score_matched_trimmed$CountTargets<numberOfTargets_Thr
Positive_Set_DrugGene_Pairs=na.omit(data.frame(Drugbank_Gene=sab_score_matched_trimmed$Targets[cond1], 
                                               Chemical.Name=sab_score_matched_trimmed$ProbeName[cond1]))
Positive_Set_DrugGene_Pairs=unique(Positive_Set_DrugGene_Pairs)
Positive_Set_DrugGene_Pairs=Positive_Set_DrugGene_Pairs[Positive_Set_DrugGene_Pairs$Drugbank_Gene %in%
                                                          rownames(onTarget$corrMat_bothScreens),]
Positive_Set_DrugGene_Pairs$Avg_Rating=sab_score_matched$`Avg Rating (in cells)`[match(Positive_Set_DrugGene_Pairs$Chemical.Name, sab_score_matched$`Probe Name`)]
Positive_Set_DrugGene_Pairs=Positive_Set_DrugGene_Pairs[Positive_Set_DrugGene_Pairs$Avg_Rating>Avg_Rating_thr,1:2]
print(dim(Positive_Set_DrugGene_Pairs))
set.seed(seedNumber)
Negative_Set_DrugGene_Pairs=randomize_DF_with_unique_Pairs_V2(Positive_Set_DrugGene_Pairs)
Positive_Set_DrugGene_Pairs$Label=1
Negative_Set_DrugGene_Pairs$Label=0
Complete_Set=rbind(Positive_Set_DrugGene_Pairs, Negative_Set_DrugGene_Pairs)
Complete_Set$Label= factor(Complete_Set$Label, labels = c('Negative', 'Positive'))
Complete_Set$drug_BroadID=onTarget$drugCategory$broad_id_trimmed[match(tolower(Complete_Set$Chemical.Name), tolower(onTarget$drugCategory$name))]
Complete_Set$crispr_corr=sapply(1:nrow(Complete_Set), function(x)
  err_handle(onTarget$corrMat_zscore[as.character(Complete_Set$Drugbank_Gene)[x], as.character(Complete_Set$drug_BroadID)[x]]))
Complete_Set$shRNA_corr=sapply(1:nrow(Complete_Set), function(x)
  err_handle(onTarget$corrMat_shRNA_zscore[as.character(Complete_Set$Drugbank_Gene)[x], as.character(Complete_Set$drug_BroadID)[x]]))
Complete_Set$both_corr=sapply(1:nrow(Complete_Set), function(x)
  err_handle(onTarget$corrMat_bothScreens_zscore[as.character(Complete_Set$Drugbank_Gene)[x], as.character(Complete_Set$drug_BroadID)[x]]))
Complete_Set=na.omit(Complete_Set)
roc_curve_crispr=roc(Complete_Set$Label, Complete_Set$crispr_corr)$auc
roc_curve_shRNA=roc(Complete_Set$Label, Complete_Set$shRNA_corr)$auc
roc_curve_both=roc(Complete_Set$Label, Complete_Set$both_corr)$auc
c(roc_curve_crispr, roc_curve_shRNA, roc_curve_both)
Complete_Set_A4=Complete_Set

# A1
Citations_Count_thr=4; seedNumber=1; Count_KnownTarget_thr=2
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
Complete_Set$crispr_corr=sapply(1:nrow(Complete_Set), function(x)
  err_handle(onTarget$corrMat_zscore[as.character(Complete_Set$Drugbank_Gene)[x], as.character(Complete_Set$drug_BroadID)[x]]))
Complete_Set$shRNA_corr=sapply(1:nrow(Complete_Set), function(x)
  err_handle(onTarget$corrMat_shRNA_zscore[as.character(Complete_Set$Drugbank_Gene)[x], as.character(Complete_Set$drug_BroadID)[x]]))
Complete_Set$both_corr=sapply(1:nrow(Complete_Set), function(x)
  err_handle(onTarget$corrMat_bothScreens_zscore[as.character(Complete_Set$Drugbank_Gene)[x], as.character(Complete_Set$drug_BroadID)[x]]))
Complete_Set=na.omit(Complete_Set)
Complete_Set_A1=Complete_Set

# TOgether
Complete_Set_both=rbind(Complete_Set_A1, Complete_Set_A4)
roc_curve_crispr=roc(Complete_Set_both$Label, Complete_Set_both$crispr_corr)$auc
roc_curve_shRNA=roc(Complete_Set_both$Label, Complete_Set_both$shRNA_corr)$auc
roc_curve_both=roc(Complete_Set_both$Label, Complete_Set_both$both_corr)$auc
c(roc_curve_crispr, roc_curve_shRNA, roc_curve_both)

Complete_Set_both$Chemical.Name

ggplot(Complete_Set_both, aes(x=Label, y=both_corr))+
  geom_boxplot()

