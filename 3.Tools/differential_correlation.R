# Perform differential correlation test of Predicted target by known target activity
require(matrixStats)
require(statar)
###########################################################################
# create drug list matrix which potentially have an alternative/alterantive target
###########################################################################
df2plot=onTarget$PredvsKnown_scores_bothScreens
df2plot$KnownTarget_corrMax_modified=df2plot$KnownTarget_corrMax
df2plot$KnownTarget_corrMax_modified[is.na(df2plot$KnownTarget_corrMax_modified)]=0
df2plot$KnownTarget_corrMax_modified[is.infinite(df2plot$KnownTarget_corrMax_modified)]=0
df2plot$Score=as.numeric(as.character(df2plot$Score))

df2plot$DrugsofInterest=df2plot$Score - df2plot$KnownTarget_corrMax_modified > 0.2 & df2plot$Score>0.3
df2plot[df2plot$DrugsofInterest,]
df2plot$DrugsofInterest_Names=''

df2plot$DrugsofInterest_Names[df2plot$DrugsofInterest]=
  paste(as.character(df2plot$CommonDrugName[df2plot$DrugsofInterest]),
        as.character(df2plot$Best_among_KnownTarget_based_onCorr[df2plot$DrugsofInterest]),
        as.character(df2plot$PredTarget[df2plot$DrugsofInterest]),
        sep=';')
###########################################################################
# Write the drugs of interest
###########################################################################
df2plot=as.data.table(df2plot)
colnames(df2plot)
drugsofInterest=df2plot[DrugsofInterest & drugCategory=='targeted cancer',c(1, 3, 10, 5)]
write.csv(drugsofInterest, '2.Data/drugsofInterest.csv')
###########################################################################
# Expression distribution of targets
###########################################################################
expression_Known_Target=do.call(rbind, lapply(drugsofInterest$Best_among_KnownTarget_based_onCorr,
                                              function(x) err_handle(data.table(Exp=onTarget$expression[as.character(x),], gene=as.character(x))))[-6])
Exp_median_Targets=aggregate(Exp ~ gene, data = expression_Known_Target, function(x) median(x))
Background_Exp=median(rowMedians(onTarget$expression))
Exp_median_Targets$drugName=drugsofInterest$CommonDrugName[match(Exp_median_Targets$gene, drugsofInterest$Best_among_KnownTarget_based_onCorr)]
Exp_median_Targets$annotation=paste(Exp_median_Targets$gene, Exp_median_Targets$drugName, sep = '\n')

tiff('4.Results/DOI_KnownTarget_Expression.tiff')
ggplot(data=Exp_median_Targets, aes(x=annotation, y=Exp)) +
  geom_bar(stat="identity")+
  geom_hline(yintercept = Background_Exp)+
  labs(x='Target', y='Median Expression in TPM')+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
dev.off()

# I conclude that only BTK (ibrutinib) and CSF1R (tandutinib) 
# genes are not expressed at all
###########################################################################
# Drugs whose Primary targets are not Expressed
###########################################################################
lowTargetExpression_Drugs=Exp_median_Targets$drugName[Exp_median_Targets$Exp < 1]
###########################################################################
# Confirming whether it is true for all the targets
###########################################################################
Targets_tandutinib=strsplit(as.character(unlist(drugsofInterest[match(as.character(
  lowTargetExpression_Drugs[2]), drugsofInterest$CommonDrugName),4])),
  ', ')[[1]]
Targets_tandutinib=c(Targets_tandutinib, 'PDGFRA')

expression_Known_Target=do.call(rbind, lapply(Targets_tandutinib,
                                              function(x) err_handle(data.table(Exp=onTarget$expression[as.character(x),], gene=as.character(x)))))
Exp_median_Targets=aggregate(Exp ~ gene, data = expression_Known_Target, function(x) median(x))
Background_Exp=median(rowMedians(onTarget$expression))
Exp_median_Targets$drugName='tandutinib'
Exp_median_Targets$annotation=paste(Exp_median_Targets$gene, Exp_median_Targets$drugName, sep = '\n')
tiff('4.Results/tandutinib_KnownvsPredTarget_Expression.tiff')
ggplot(data=Exp_median_Targets, aes(x=annotation, y=Exp)) +
  geom_bar(stat="identity")+
  geom_hline(yintercept = Background_Exp)+
  geom_hline(yintercept = 1, color='red', type='dotted')+
  labs(x='Target', y='Median Expression in TPM')+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
dev.off()
###########################################################################
# Identifying drugs 
###########################################################################
test_differential_correlation<-function(infunc_drugName='ibrutinib',
                                        infunc_PredTarget='EGFR',
                                        infunc_KnownTarget='BTK',
                                        median_exp=Exp_median_Targets$Exp[Exp_median_Targets$drugName=='ibrutinib']){
  drugID=onTarget$drugCategory$broad_id_trimmed[match(infunc_drugName,
                                                      onTarget$drugCategory$name)]
  common_Columns=Reduce(intersect, list(a=colnames(onTarget$avana),
                                        b=colnames(onTarget$drug_prism),
                                        c=colnames(onTarget$expression)))
  if(median_exp<1){
    infunc_KnownTarget_highID=common_Columns[onTarget$expression[infunc_KnownTarget, common_Columns]>1]
    infunc_KnownTarget_lowID=common_Columns[onTarget$expression[infunc_KnownTarget, common_Columns]<1]
  } else{
    infunc_KnownTarget_tiled=xtile(onTarget$expression[infunc_KnownTarget, common_Columns], 3)
    infunc_KnownTarget_lowID=common_Columns[infunc_KnownTarget_tiled==1]
    infunc_KnownTarget_highID=common_Columns[infunc_KnownTarget_tiled==3]
  }
  common_Columns_avanaPrism=Reduce(intersect, list(a=colnames(onTarget$avana),
                                        b=colnames(onTarget$drug_prism)))
  
  rbind(PredCorr=cor.test(onTarget$avana[infunc_PredTarget,common_Columns_avanaPrism],
                          onTarget$drug_prism[drugID,common_Columns_avanaPrism], method = 's')[c(3,4)],
        
        PredCorr_lowKnown=cor.test(onTarget$avana[infunc_PredTarget,infunc_KnownTarget_lowID],
                                   onTarget$drug_prism[drugID,infunc_KnownTarget_lowID], method = 's')[c(3,4)],
        
        PredCorr_HighKnown=cor.test(onTarget$avana[infunc_PredTarget,infunc_KnownTarget_highID],
                                    onTarget$drug_prism[drugID,infunc_KnownTarget_highID], method = 's')[c(3,4)],
        
        KnownCorr=cor.test(onTarget$avana[infunc_KnownTarget,common_Columns_avanaPrism],
                           onTarget$drug_prism[drugID,common_Columns_avanaPrism], method = 's')[c(3,4)],
        
        KnownCorr_lowKnown=cor.test(onTarget$avana[infunc_KnownTarget,infunc_KnownTarget_lowID],
                                    onTarget$drug_prism[drugID,infunc_KnownTarget_lowID], method = 's')[c(3,4)],
        
        KnownCorr_HighKnown=cor.test(onTarget$avana[infunc_KnownTarget,infunc_KnownTarget_highID],
                                     onTarget$drug_prism[drugID,infunc_KnownTarget_highID], method = 's')[c(3,4)]
  )
  
}
drugsofInterest=drugsofInterest[,-5]
diff_corrv2=apply(drugsofInterest, 1, function(x)
  err_handle(test_differential_correlation(as.character(x[1]),
                                           as.character(x[2]),
                                           as.character(x[3]),
                                           Exp_median_Targets$Exp[Exp_median_Targets$gene==
                                                                    as.character(x[3])])) )

names(diff_corrv2)=as.character(drugsofInterest$CommonDrugName)
diff_corr; K=17
drugsofInterest[K,]
diff_corr[[K]]
################################################################################
# If Drug response variance and Target response variance is confouding our results
################################################################################
variance_DrugResponse=apply(onTarget$drug_prism, 1, function(x) var(na.omit(x)))
background_Variance_drugResp=median(variance_DrugResponse)
drugsofInterest$variance_DrugResponse=variance_DrugResponse[
  match(drugsofInterest$CommonDrugName, onTarget$drugCategory$name)]
variance_KOResponse=apply(onTarget$avana, 1, function(x) var(na.omit(x)))
background_variance_KOResponse=median(variance_KOResponse)
drugsofInterest$variance_KOResponse=variance_KOResponse[
  match(drugsofInterest$Best_among_KnownTarget_based_onCorr, names(variance_KOResponse))]
drugsofInterest[drugsofInterest$variance_KOResponse<background_variance_KOResponse,]
drugsofInterest[drugsofInterest$variance_DrugResponse<background_Variance_drugResp,]
