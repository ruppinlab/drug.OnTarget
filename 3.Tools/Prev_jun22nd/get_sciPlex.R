# match drug response i have
###########################################################################
# Reading Single Cell DrugPerturb-seq
###########################################################################
require(data.table)
require(parallel)
setwd('/Users/sinhas8/Project_OffTarget/2.Data/')
# PossibleRowName=sciPlex[,1]
sciPlex=fread('Supplementary_Table_5.txt', skip=1, fill=T, sep=',')
ace2_ID=grep('TMPRSS2',
             sciPlex$`A549	S1133	ENSG00000000003.14	TSPAN6	OK	(Intercept)	-3.5775032	0.172417	-20.7491	2.67e-88	0	count	3.7600275e-83	72hours`)
sciPlex_ace2=sciPlex[ace2_ID,]
dim(sciPlex_ace2)
sciPlex_ace2_list=apply(sciPlex_ace2, 1, function(x) strsplit(x, '\t'))
names(sciPlex_ace2_list[[1]])<-NULL
sciPlex_ace2_list_len=unlist(mclapply(sciPlex_ace2_list, function(x) length(x[[1]]), mc.cores = 4))
overlapping_only_id=unlist(mclapply(sciPlex_ace2_list, function(x) 
  x[[1]][2] %in% overlapping_drugs_id, mc.cores=4))
ace2_scRNA=sapply(sciPlex_ace2_list[overlapping_only_id], function(x) x[[1]])
ace2_scRNA[25:30]
sapply(ace2_scRNA, function(x) x[4]=='ACE2')

setwd('/Users/sinhas8/Project_OffTarget/2.Data/')
meta_sciPlex=read.csv('/Users/sinhas8/Downloads/aax6234-Srivatsan-Table-S3.txt', sep='\t')
onTarget=readRDS('/Users/sinhas8/Project_OffTarget/2.Data/onTarget.RDS')
overlapping_drugs_id=levels(meta_sciPlex$catalog_number)[
  !is.na(match(tolower(levels(meta_sciPlex$name)),
               tolower(onTarget$drugsCommonName)))]
overlapping_drugs=levels(meta_sciPlex$name)[!is.na(match(tolower(levels(meta_sciPlex$name)),
                                                                      tolower(onTarget$drugsCommonName)))]

head(meta_sciPlex)
overlapping_drugs