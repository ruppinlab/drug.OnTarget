# Validation A9 - COSMICv91 Drug Resistance file
require(data.table)
require(pROC)
setwd('/Users/sinhas8/Project_OffTarget/')
source('drug.OnTarget/3.Tools/myCustom_functions.R')
CosmicResistanceMutations=read.table('CosmicResistanceMutations.tsv', sep='\t')
colnames(CosmicResistanceMutations)=as.character(unlist(CosmicResistanceMutations[1,]))
head(sort(table(CosmicResistanceMutations$`Pubmed Id`), decreasing = T))

head(CosmicResistanceMutations)
CosmicResistanceMutations=CosmicResistanceMutations[-1,]
CosmicResistanceMutations=CosmicResistanceMutations[-1,]
table(CosmicResistanceMutations$`Drug Name`)
CosmicResistanceMutations[CosmicResistanceMutations$`Drug Name`=='Ibrutinib',]
CosmicResistanceMutations$trimmed_geneName= sapply(CosmicResistanceMutations$`Gene Name`,function(x) strsplit(as.character(x), '_')[[1]][1])
DrugTargetPair= data.frame(drugname=CosmicResistanceMutations$`Drug Name`, GeneName=CosmicResistanceMutations$trimmed_geneName)
Evidence_Strength=table(paste(DrugTargetPair$drugname, DrugTargetPair$GeneName))
# Evidence_Strength=table(Evidence_Strength)
DrugTargetPair=unique(DrugTargetPair)
DrugTargetPair$drugName_geneName=paste(DrugTargetPair$drugname, DrugTargetPair$GeneName)

DrugTargetPair=cbind(DrugTargetPair,
                     Evidence_Strength=as.numeric(Evidence_Strength[
                       match(DrugTargetPair$drugName_geneName, names(Evidence_Strength))]))
numberofTargets=table(DrugTargetPair$drugname)
DrugTargetPair$numberofTargets=numberofTargets[match(DrugTargetPair$drugname, names(numberofTargets))]

test_AUC<-function(seedNumber=1, numberOfTargets_Thr=Inf, Evidence_Strength_thr=0){
  cond1=DrugTargetPair$numberofTargets<numberOfTargets_Thr
  cond2=DrugTargetPair$Evidence_Strength> Evidence_Strength_thr
  Positive_Set_DrugGene_Pairs=na.omit(data.frame(Drugbank_Gene=as.character(DrugTargetPair$GeneName)[cond1 & cond2], 
                                                 Chemical.Name=as.character(DrugTargetPair$drugname)[cond1 & cond2]))
  
  Positive_Set_DrugGene_Pairs=unique(Positive_Set_DrugGene_Pairs)
  dim(Positive_Set_DrugGene_Pairs)
  Positive_Set_DrugGene_Pairs=Positive_Set_DrugGene_Pairs[Positive_Set_DrugGene_Pairs$Drugbank_Gene %in%
                                                            rownames(onTarget$corrMat_bothScreens),]
  dim(Positive_Set_DrugGene_Pairs)
  Positive_Set_DrugGene_Pairs=Positive_Set_DrugGene_Pairs[tolower(Positive_Set_DrugGene_Pairs$Chemical.Name)  %in%
                                                            tolower(onTarget$drugCategory_prism2$name),]
  print(dim(Positive_Set_DrugGene_Pairs))
  # print(length(table(as.character(Positive_Set_DrugGene_Pairs$Chemical.Name))))
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
  Complete_Set$prism2=sapply(1:nrow(Complete_Set), function(x)
    err_handle(onTarget$corrMat_prism2[as.character(Complete_Set$Drugbank_Gene)[x], as.character(Complete_Set$drug_BroadID)[x]]))
  Complete_Set=na.omit(Complete_Set)
  
  roc_curve_crispr=roc(Complete_Set$Label, Complete_Set$crispr_corr)$auc
  roc_curve_shRNA=roc(Complete_Set$Label, Complete_Set$shRNA_corr)$auc
  roc_curve_both=roc(Complete_Set$Label, Complete_Set$both_corr)$auc
  roc_curve_prism2=roc(Complete_Set$Label, Complete_Set$prism2)$auc
  c(roc_curve_crispr, roc_curve_shRNA, roc_curve_both, roc_curve_prism2)
  # roc_obj=roc(Complete_Set$Label, Complete_Set$both_corr)
  # range(Complete_Set$both_corr[Complete_Set$Label=='Positive'])
  # ggroc(roc_obj)+
  #   theme_minimal() + ggtitle("My ROC curve") + 
  #   geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed")
  
}
sapply(0, function(x) rowMeans(sapply(1:10, function(y) test_AUC(seedNumber=y, numberOfTargets_Thr=Inf, Evidence_Strength_thr=x))) )


# Drugs whose annotation is wrong
# match(tolower(as.character(df2plot$CommonDrugName[df2plot$DrugsofInterest])), tolower(DrugTargetPair$drugname))
