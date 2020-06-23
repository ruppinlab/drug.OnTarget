# A6 Validation including hits from Sheltzer natCan 2020 and STM 2019.
require(ggrepel)
# Target Identification of OTS964

sum(onTarget$corrMat_bothScreens[,match('OTS167', onTarget$drugCategory$name)]>0.19)
# Not available in our Screens

# All drugs used in the study
sheltzer_STM_hits=readxl::read_xlsx('2.Data/Sheltzer_STM.xlsx')
sheltzer_STM_hits_matched=onTarget$drugCategory$name[na.omit(match(gsub('-','',tolower(sheltzer_STM_hits$Drug)), 
                                    gsub('-','',tolower(onTarget$drugCategory$name))))]

onTarget$PredvsKnown_scores[match(sheltzer_STM_hits_matched,
                                  onTarget$PredvsKnown_scores$CommonDrugName),]

df2plot=onTarget$PredvsKnown_scores
df2plot$DrugsofInterest=df2plot$KnownTarget_corrMax < 0.1 & df2plot$Score > 0.3
df2plot$DrugsofInterest_Names=''
df2plot$DrugsofInterest_Names[df2plot$DrugsofInterest]=as.character(df2plot$CommonDrugName[df2plot$DrugsofInterest])

df2plot$DrugsofInterest_Names[df2plot$DrugsofInterest]=
  paste(as.character(df2plot$CommonDrugName[df2plot$DrugsofInterest]),
        as.character(df2plot$Best_among_KnownTarget_based_onCorr[df2plot$DrugsofInterest]),
        as.character(df2plot$PredTarget[df2plot$DrugsofInterest]),
        sep=';')


ggplot(df2plot, aes(x=KnownTarget_corrMax, y=Score, color=DrugsofInterest,
                                        label=DrugsofInterest_Names))+
  geom_point()+
  geom_abline(intercept = 0, slope = 1)+
  geom_label_repel(nudge_y = 0.25, nudge_x = 0.15)