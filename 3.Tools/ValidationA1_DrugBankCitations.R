# ValidataionA1 - By Curating BioGrid Chemical Interaction database we would like to test 
require(stats)
require(reshape2)
require(statar)
require(ggplot2)
biogrid=read.csv('/Users/sinhas8/Project_OffTarget/2.Data/BIOGRID-CHEMICALS-3.5.182.chemtab.txt', sep='\t')
common_Drug_lowerCase=tolower(levels(biogrid$Chemical.Name))[
  !is.na(match(tolower(levels(biogrid$Chemical.Name)), tolower(onTarget$PredvsKnown_scores$CommonDrugName) ))]
biogrid_matched=biogrid[!is.na(match(tolower(biogrid$Chemical.Name), common_Drug_lowerCase)),]
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

DrugBankScore_vsOurScore=data.frame(biogrid_matched_Subsetlist_withScore, 
                                    PredvsKnown_scores_matched[match(
                                      tolower(biogrid_matched_Subsetlist_withScore$Chemical.Name),
                                      tolower(PredvsKnown_scores_matched$CommonDrugName)),])

# Only testing the hypothesis on the targeted therapy drugs (relevant subset)
targeted_drugs=onTarget$drugCategory$broad_id[onTarget$drugCategory$drug_category=='targeted cancer']
targeted_drugs_trimmed=sapply(as.character(targeted_drugs), function(x) strsplit(x, '-')[[1]][2])
DrugBankScore_vsOurScore_CancerDrugs=DrugBankScore_vsOurScore[!is.na(match(DrugBankScore_vsOurScore$drugName,
                                                                           targeted_drugs_trimmed)),]
DrugBankScore_vsOurScore_CancerDrugs[which(DrugBankScore_vsOurScore_CancerDrugs$KnownTarget_corrRank_min==1),]
DrugBankScore_vsOurScore_CancerDrugs$Count_KnownTarget=sapply(DrugBankScore_vsOurScore_CancerDrugs$KnownTarget,
                                                              function(x) length(strsplit(as.character(x), '\\,')[[1]])) 
singleTargetDrugs=which(DrugBankScore_vsOurScore_CancerDrugs$Count_KnownTarget<4 )
inhibitors=which(DrugBankScore_vsOurScore_CancerDrugs$Action=='inhibitor')
ggplot(DrugBankScore_vsOurScore_CancerDrugs[intersect(singleTargetDrugs, inhibitors),], 
       aes(x= factor(xtile(Score, 3)),
           y= KnownTarget_corrRank_min))+
  geom_boxplot()
