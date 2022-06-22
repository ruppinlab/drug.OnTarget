drugList_ofInterest=onTarget$drugCategory[which(onTarget$drugCategory$drug_category == 'targeted cancer'),]
write.table(drugList_ofInterest,
            '/Users/sinhas8/Project_OffTarget/2.Data/drugList_ofInterest.tsv',
            sep='\t')

# For Alejandro
load('/Users/sinhas8/Downloads/drugs2targets.RData')
