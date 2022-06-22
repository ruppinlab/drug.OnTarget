###########################################################################
# matched drug response from a technique called DrugPerturb-seq
###########################################################################
setwd('/Users/sinhas8/Project_OffTarget/2.Data/')
depmap=get(load('depmap_19Q3.RData'))
meta_sciPlex=read.csv('/Users/sinhas8/Downloads/aax6234-Srivatsan-Table-S3.txt', sep='\t')
overlapping_drugs_id=levels(meta_sciPlex$catalog_number)[!is.na(match(tolower(levels(meta_sciPlex$name)),
                                                                      tolower(depmap_19Q3$drug_target_map$name)))]
overlapping_drugs=levels(meta_sciPlex$name)[!is.na(match(tolower(levels(meta_sciPlex$name)),
                                                         tolower(depmap_19Q3$drug_target_map$name)))]
###########################################################################
# Reading Single Cell DrugPerturb-seq
###########################################################################
require(data.table)
setwd('/Users/sinhas8/Project_OffTarget/2.Data/')
# PossibleRowName=sciPlex[,1]
sciPlex=fread('Supplementary_Table_5.txt', skip=1, fill=T, sep=',')
sciPlex_list=apply(sciPlex, 1, function(x) strsplit(x, '\t'))
head(names(sciPlex_list[[1]]))
names(sciPlex_list[[1]])<-NULL
sciPlex_list_len=unlist(mclapply(sciPlex_list, function(x) length(x[[1]]), mc.cores = 4))
overlapping_only_id=unlist(mclapply(sciPlex_list, function(x) 
  x[[1]][2] %in% overlapping_drugs_id, mc.cores=4))
sciPlex_list_ovp=sciPlex_list[overlapping_only_id]
table(sapply(sciPlex_list_ovp, length))
###########################################################################
# 
###########################################################################
sciPlex_test1=fread('Supplementary_Table_5.txt', skip=1, fill=T, sep='\t', nrows = 1120480)
sciPlex_test2=fread('Supplementary_Table_5.txt', skip=1120481, fill=T, sep='\t', nrows=1120491)
tail(sciPlex_test2, 10)