# A6 Validation including hits from Sheltzer natCan 2020 and STM 2019.
require(ggrepel)
# All drugs used in the study
sheltzer_STM_hits=readxl::read_xlsx('2.Data/Sheltzer_STM.xlsx')
sheltzer_STM_hits_matched=onTarget$drugCategory$name[na.omit(match(gsub('-','',tolower(sheltzer_STM_hits$Drug)), 
                                                                   gsub('-','',tolower(onTarget$drugCategory$name))))]
onTarget$PredvsKnown_scores[na.omit(match(gsub('-','',tolower(sheltzer_STM_hits$Drug)),
                                         gsub('-','',tolower(onTarget$PredvsKnown_scores$CommonDrugName)))),]
###########################################################################
# OTS167
###########################################################################
sum(onTarget$corrMat_bothScreens[,match('OTS167', onTarget$drugCategory$name)]>0.19)
# Not available in our Screens
###########################################################################
# AZ3146
###########################################################################
AZ3146_score=onTarget$corrMat_bothScreens[,grep('AZ3146',
                                                onTarget$PredvsKnown_scores_bothScreens$CommonDrugName, ignore.case = T)]
AZ3146_score=sort(AZ3146_score, decreasing = T)
AZ3146_score[grep('TTK',names(AZ3146_score))]
# Our pipeline suggests a diff target in comparison to Sheltzer et al 
###########################################################################
# nutlin
###########################################################################
nutlin_score=onTarget$corrMat_bothScreens[,grep('nutlin-3',
                                                onTarget$PredvsKnown_scores_bothScreens$CommonDrugName, ignore.case = T)]
nutlin_score=sort(nutlin_score, decreasing = T)
grep('MDM2',names(nutlin_score))
grep('TP53',names(nutlin_score))
###########################################################################
# rapamycin
###########################################################################
rapamycin_score=onTarget$corrMat_bothScreens[,grep('rapamycin',
                                                onTarget$PredvsKnown_scores_bothScreens$CommonDrugName, ignore.case = T)]
rapamycin_score=sort(rapamycin_score, decreasing = T)
grep('MDM2',names(rapamycin_score))
grep('TP53',names(rapamycin_score))

###########################################################################
# OTS964
###########################################################################
sum(onTarget$corrMat_bothScreens[,match('OTS964', onTarget$drugCategory$name)]>0.19)
# Not available in our Screens
###########################################################################
# For the 10 drugs used in the study:: Can predict for four
###########################################################################
df2plot=onTarget$PredvsKnown_scores_bothScreens
df2plot$KnownTarget_corrMax_modified=df2plot$KnownTarget_corrMax
df2plot$KnownTarget_corrMax_modified[is.na(df2plot$KnownTarget_corrMax_modified)]=0
df2plot$KnownTarget_corrMax_modified[is.infinite(df2plot$KnownTarget_corrMax_modified)]=0
df2plot$Score=as.numeric(as.character(df2plot$Score))

df2plot$DrugsofInterest=df2plot$Score - df2plot$KnownTarget_corrMax_modified > 0.2 & df2plot$Score>0.3

df2plot$DrugsofInterest_Names=''

df2plot$DrugsofInterest_Names[df2plot$DrugsofInterest]=
  paste(as.character(df2plot$CommonDrugName[df2plot$DrugsofInterest]),
        as.character(df2plot$Best_among_KnownTarget_based_onCorr[df2plot$DrugsofInterest]),
        as.character(df2plot$PredTarget[df2plot$DrugsofInterest]),
        sep=';')
tiff('4.Results/KnownvsBest.tiff', width=1200, height=800)
ggplot(df2plot, aes(x=KnownTarget_corrMax, y=Score, color=DrugsofInterest,
                                        label=DrugsofInterest_Names))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1)+
  geom_label_repel(nudge_y = 0.25, nudge_x = 0.15)+
  theme_bw(base_size = 25)+
  labs(x='Corr Rho of Known Target', y='Rho of Best Hit')+
  theme(legend.position = 'none')+
  ylim(c(-0.15, 0.8))+
  xlim(c(-0.15, 0.8))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)
dev.off()

tiff('4.Results/KnownvsBest_Targeted_Therapy.tiff', width=800, height=800)
ggplot(df2plot[df2plot$drugCategory=='targeted cancer',], aes(x=KnownTarget_corrMax, y=Score, color=DrugsofInterest,
                    label=DrugsofInterest_Names))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1)+
  geom_label_repel(nudge_y = 0.25, nudge_x = 0.15)+
  theme_bw(base_size = 25)+
  labs(x='Corr Rho of Known Target', y='Rho of Best Hit')+
  theme(legend.position = 'none')+
  # ylim(c(-0.15, 0.8))+
  # xlim(c(-0.15, 0.8))+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)
dev.off()


###########################################################################
# From Sheltzer NatCan 2020
###########################################################################
# BI-2536
onTarget$PredvsKnown_scores_bothScreens[grep('BI-2536',onTarget$PredvsKnown_scores_bothScreens$CommonDrugName, ignore.case = T),]
# bortezomib
onTarget$PredvsKnown_scores_bothScreens[grep('bortezomib',onTarget$PredvsKnown_scores_bothScreens$CommonDrugName, ignore.case = T),]
# cysteine
onTarget$PredvsKnown_scores_bothScreens[grep('cysteine',onTarget$PredvsKnown_scores_bothScreens$CommonDrugName, ignore.case = T),]

# Find approved drugs among top hits
DOI=df2plot$CommonDrugName[df2plot$DrugsofInterest & df2plot$drugCategory=='targeted cancer']
FDAapproved_drugs$name[unlist(sapply(DOI, function(x) grep(x, FDAapproved_drugs$name, ignore.case = T)))]

