onTarget_comp=readRDS('/Users/sinhas8/Project_OffTarget/2.Data/onTarget_v7.RDS')
id_matched_lin=c(grep('1776',onTarget_comp$drugsCommonName),
                 grep('OTS964',onTarget_comp$drugsCommonName, ignore.case = T),
                 grep('OTS514',onTarget_comp$drugsCommonName, ignore.case = T),
                 
                 grep('03758309',onTarget_comp$drugsCommonName),
                 # 
                 grep('SCIO469',onTarget_comp$drugsCommonName, ignore.case = T),
                 # Ralimetinib
                 grep('2228820',onTarget_comp$drugsCommonName, ignore.case = T),
                 # Ricolinostat
                 grep('ACY-1215',onTarget_comp$drugsCommonName, ignore.case = T),
                 grep('Citarinostat',onTarget_comp$drugsCommonName, ignore.case = T),
                 grep('PAC-1',onTarget_comp$drugsCommonName, ignore.case = T),
                 grep('K30064966',colnames(onTarget_comp$corrMat), ignore.case = T)
                 )
Lin_comparison=list(cor_profile=onTarget_comp$corrMat[,id_matched_lin],
                    drugName=onTarget_comp$drugsCommonName[id_matched_lin],
                    knownTarget=onTarget_comp$Annotated_Target[id_matched_lin],
                    PredTarget=onTarget_comp$Top_predicted_Target_bothScreens[id_matched_lin,],
                    toPlot=onTarget_comp$PredvsKnown_scores[colnames(onTarget_comp$corrMat)[id_matched_lin],]
                    )



df2export=Lin_comparison$toPlot

write.csv(df2export,
          '/Users/sinhas8/Project_OffTarget/4.Results/comparison2Lin.csv')

onTarget_comp$drugsCommonName[id_matched_lin]

df2plot=Lin_comparison$toPlot

data.frame(df2plot[1:3,c(1, 2, 3)], type='Pred')
