# Validation in another dataset - From Miquel dataset
source('drug.OnTarget/3.Tools/myCustom_functions.R')
require(data.table)
gs2=read.csv('/Users/sinhas8/Downloads/drugbank_chembldrug.tsv', sep='\t')
gs2=gs2[gs2$direction == -1,]
inchimapping=readxl::read_xlsx('/Users/sinhas8/Downloads/DSSTox_Identifiers_and_CASRN.xlsx')
mapping=match(gs2$inchikey, inchimapping$standard_InChIKey)
names_mapped=inchimapping[mapping,]
inchimapping<-NULL
gs2$CommonName=names_mapped$preferred_name
gs2$CommonName=tolower(gs2$CommonName)
gs2_matched=gs2[!is.na(match(tolower(gs2$CommonName), tolower(onTarget$drugsCommonName))),]

uniProt_geneName_map=read.csv('/Users/sinhas8/Downloads/geneName_to_UniprotSwissID.txt', sep='\t')
gs2_matched$geneName=uniProt_geneName_map$HGNC.symbol[match(gs2_matched$uniprot_ac,uniProt_geneName_map$UniProtKB.Swiss.Prot.ID)]
Count_Targets=table(gs2_matched$CommonName)
gs2_matched$Count_Targets_allType=Count_Targets[match(gs2_matched$CommonName, names(Count_Targets))]
gs2_matched=na.omit(gs2_matched)
Count_Targets=table(gs2_matched$CommonName)
gs2_matched$Count_Targets_ProteinOnly=Count_Targets[match(gs2_matched$CommonName, names(Count_Targets))]

gs2_matched_trimmed=data.frame(ProbeName=gs2_matched$CommonName,
                                     Targets=gs2_matched$geneName)
gs2_matched_trimmed$Targets=as.character(gs2_matched_trimmed$Targets)
gs2_matched_trimmed$CountTargets=gs2_matched$Count_Targets_ProteinOnly


test_AUC<-function(seedNumber=1, numberOfTargets_Thr=2){
  cond1=gs2_matched_trimmed$CountTargets<numberOfTargets_Thr
  Positive_Set_DrugGene_Pairs=na.omit(data.frame(Drugbank_Gene=gs2_matched_trimmed$Targets[cond1], 
                                                 Chemical.Name=gs2_matched_trimmed$ProbeName[cond1]))
  Positive_Set_DrugGene_Pairs=unique(Positive_Set_DrugGene_Pairs)
  Positive_Set_DrugGene_Pairs=Positive_Set_DrugGene_Pairs[Positive_Set_DrugGene_Pairs$Drugbank_Gene %in%
                                                            rownames(onTarget$corrMat_bothScreens),]
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
  roc_curve_crispr=roc(Complete_Set$Label, Complete_Set$crispr_corr)$auc
  roc_curve_shRNA=roc(Complete_Set$Label, Complete_Set$shRNA_corr)$auc
  roc_curve_both=roc(Complete_Set$Label, Complete_Set$both_corr)$auc
  c(roc_curve_crispr, roc_curve_shRNA, roc_curve_both)
  # roc_obj=roc(Complete_Set$Label, Complete_Set$both_corr)
  # range(Complete_Set$both_corr[Complete_Set$Label=='Positive'])
  # ggroc(roc_obj)+
  #   theme_minimal() + ggtitle("My ROC curve") + 
  #   geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed")
  
}
test_AUC(seedNumber=1, numberOfTargets_Thr=2)
