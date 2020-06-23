# A5 Validation using ChemicalProbes.org Scientific Advisory Board recommendation score
setwd('/Users/sinhas8/Project_OffTarget/')
sab_score=readxl::read_xlsx('2.Data/SAB_Score_chemical_Probes.xlsx')
sab_score$`Probe Name`=tolower(sab_score$`Probe Name`)
sab_score_matched=sab_score[!is.na(match(tolower(sab_score$`Probe Name`), tolower(onTarget$drugsCommonName))),]
sab_score_matched$CountTargets=sapply(sab_score_matched$`Protein target`,
                                      function(x) length(strsplit(x, ', ')[[1]]))
sab_score_matched_trimmed=data.frame(ProbeName=rep(sab_score_matched$`Probe Name`, sab_score_matched$CountTargets),
                                     Targets=unlist(sapply(sab_score_matched$`Protein target`, function(x) strsplit(x, ', ')[[1]])))
sab_score_matched_trimmed$Targets=as.character(sab_score_matched_trimmed$Targets)
sab_score_matched_trimmed$Targets[grep('IDH1',sab_score_matched_trimmed$Targets)]='IDH1'
sab_score_matched_trimmed$Targets[grep('IDH2',sab_score_matched_trimmed$Targets)]='IDH2'
sab_score_matched_trimmed$Targets[grep('BRAF',sab_score_matched_trimmed$Targets)]='BRAF'
CountTargets=table(sab_score_matched_trimmed$ProbeName)
sab_score_matched_trimmed$CountTargets=CountTargets[match(sab_score_matched_trimmed$ProbeName,
                                                          names(CountTargets))]

test_AUC<-function(seedNumber=1, numberOfTargets_Thr=2, Avg_Rating_thr=3){
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
    err_handle(onTarget$corrMat[as.character(Complete_Set$Drugbank_Gene)[x], as.character(Complete_Set$drug_BroadID)[x]]))
  Complete_Set$shRNA_corr=sapply(1:nrow(Complete_Set), function(x)
    err_handle(onTarget$corrMat_shRNA[as.character(Complete_Set$Drugbank_Gene)[x], as.character(Complete_Set$drug_BroadID)[x]]))
  Complete_Set$both_corr=sapply(1:nrow(Complete_Set), function(x)
    err_handle(onTarget$corrMat_bothScreens[as.character(Complete_Set$Drugbank_Gene)[x], as.character(Complete_Set$drug_BroadID)[x]]))
  Complete_Set=na.omit(Complete_Set)
  ggplot(Complete_Set, aes(x=Label, y=both_corr))+
    geom_boxplot()
  roc_curve_crispr=roc(Complete_Set$Label, Complete_Set$crispr_corr)$auc
  roc_curve_shRNA=roc(Complete_Set$Label, Complete_Set$shRNA_corr)$auc
  roc_curve_both=roc(Complete_Set$Label, Complete_Set$both_corr)$auc
  c(roc_curve_crispr, roc_curve_shRNA, roc_curve_both)
}

rowMeans(sapply(1:10, function(x)
  test_AUC(seedNumber=x, numberOfTargets_Thr=2,Avg_Rating_thr=3)))


