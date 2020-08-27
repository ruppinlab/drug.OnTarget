################################################################
# Provide alejandro a set of drugs with less than four targets
################################################################
onTarget=readRDS('/Users/sinhas8/Project_OffTarget/2.Data/onTarget_v4.RDS')
drugs2targets= get(load('/Users/sinhas8/Downloads/drugs2targets.RData'))
onTarget$drugCategory$numTargets=sapply(strsplit(as.character(onTarget$drugCategory$target),', '), length)
singleTarget_Drugs=onTarget$drugCategory[onTarget$drugCategory$numTargets<2 & onTarget$drugCategory$drug_category=='targeted cancer',]
singleTarget_Drugs=singleTarget_Drugs[!is.na(singleTarget_Drugs$target),]
write.csv(singleTarget_Drugs[,-(1:2)], '/Users/sinhas8/Project_OffTarget/2.Data/singleTarget_Drugs.csv')
write.csv(onTarget$drugCategory[,-(1:2)], '/Users/sinhas8/Project_OffTarget/2.Data/All_Drugs_Targets.csv')

drug.atc= get(load('/Users/sinhas8/Downloads/drug.atc.RData'))

