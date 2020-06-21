# Drugs for Fio
fio_GeneList=read.csv('/Users/sinhas8/Project_OffTarget/2.Data/Fio_genesOfInterest.txt')
Drugs_for_a_gene<-function(GeneName_infunc){
  onTarget$Top_predicted_Drug[match(GeneName_infunc, onTarget$Top_predicted_Drug$GeneName),]
}
approved_Drugs_for_a_gene<-function(GeneName_infunc='P4HA1'){
  onTarget$Top_predicted_Drug[match(GeneName_infunc,
                                    onTarget$Top_predicted_Drug$GeneName),]
}
drugs_for_fio_GeneList=na.omit(do.call(rbind, lapply(as.character(fio_GeneList[,1]), Drugs_for_a_gene)))
drugs_for_fio_GeneList=drugs_for_fio_GeneList[order(drugs_for_fio_GeneList$Score, decreasing = T),]
drugs_for_fio_GeneList$drugsCommonName=onTarget$drugsCommonName[match(drugs_for_fio_GeneList$PredTarget,
                                                                      colnames(onTarget$corrMat))]
drugs_for_fio_GeneList$known_AnnotatedTarget=onTarget$Annotated_Target[match(drugs_for_fio_GeneList$PredTarget,
                                                                      colnames(onTarget$corrMat))]
colnames(drugs_for_fio_GeneList)[3]='ConfidenceScore'
write.table(drugs_for_fio_GeneList,
          '/Users/sinhas8/Project_OffTarget/2.Data/drugs_for_fio_GeneList.tsv',
          sep='\t', 
          row.names = F)

Drugs_for_a_gene('ACE2')

###########################################################################
# Drugs for P4HA1
###########################################################################
drugsfor_P4HA1=data.frame(Drug_BroadID=names(onTarget$corrMat['P4HA1',order(onTarget$corrMat['P4HA1',], decreasing = T)])
                          ,Correlation_Score=onTarget$corrMat['P4HA1',order(onTarget$corrMat['P4HA1',], decreasing = T)])

