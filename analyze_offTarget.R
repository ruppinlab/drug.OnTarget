###########################################################################
# Using correlation matrix, we aim to identify new targets!!
###########################################################################
require(parallel); require(ggplot2)
source('/Users/sinhas8/myCustom_functions.R')
# Files
setwd('/Users/sinhas8/Project_OffTarget/2.Data/')
depmap=get(load('depmap_19Q3.RData'))
###########################################################################
# Step 1B: Cluster and visualize the matrix
###########################################################################
cmat=readRDS('/Users/sinhas8/Project_OffTarget/2.Data/corrMat.RDS')
cmat=unlist(cmat, recursive = F)
cmat=do.call(rbind, cmat)
colnames(cmat)=rownames(prism)
rownames(cmat)=rownames(crispr)
###########################################################################
# Step 1A: Cluster and visualize the matrix
###########################################################################
plotHist<-function(testVector){
  cmat_colMaximum=data.frame(value=testVector)
  ggplot(cmat_colMaximum, aes(x=value))+
    geom_freqpoly()+
    geom_histogram(color="black", fill="white", bins=50)
}
cmat_colMax=apply(cmat, 2, max)
cmat_colMedian=apply(cmat, 2, median)
cmat_colMean=apply(cmat, 2, mean)
cmat_colSD=apply(cmat, 2, sd)
plotHist(cmat_colSD)
###########################################################################
# Step 1B: Target Correlation
###########################################################################
Targets_list=lapply(as.character(depmap$drug_target_map$target),
                    function(x) gsub(' ','',strsplit(x, '\\,')[[1]]) )
Targets_list=Targets_list[1:80]
x=80
Annotated_Target_corr=sapply(1:length(Targets_list), function(x)
  sapply(Targets_list[[x]], function(y)  err_handle(cmat[y,x]) ))
Annotated_Target_corr_mean=sapply(Annotated_Target_corr, function(x) mean(x, na.rm=T))
plotHist(Annotated_Target_corr_mean)
Annotated_Target_corr_max=sapply(Annotated_Target_corr, function(x) max(x, na.rm=T))
plotHist(Annotated_Target_corr_max)
###########################################################################
# Step 1C: Target Corr Rank 
###########################################################################
cmat_rank=apply(cmat, 2, function(x) rank(x) )
Annotated_Target_corr_rank=sapply(1:length(Targets_list), function(x)
  sapply(Targets_list[[x]], function(y)  err_handle(cmat_rank[y,x]) ))
Annotated_Target_corr_rank_mean=sapply(Annotated_Target_corr_rank, function(x) mean(x, na.rm=T))
plotHist(Annotated_Target_corr_rank_mean)
Annotated_Target_corr_rank_max=sapply(Annotated_Target_corr_rank, function(x) max(x, na.rm=T))
plotHist(Annotated_Target_corr_rank_max)
###########################################################################
# Subset of Cancer Drugs
###########################################################################
drug_details=read.csv('/Users/sinhas8/Project_OffTarget/biomarkers.csv')
drug_details_cancer=drug_details[drug_details$drug_category=='chemo' | drug_details$drug_category=='targeted cancer',]
drug_details_cancer$broad_id_trimmed=sapply(as.character(drug_details_cancer$broad_id), function(x) strsplit(x, '-')[[1]][2])
cmat_cancer=cmat[,match(drug_details_cancer$broad_id_trimmed, colnames(cmat))]
cmat_rank_cancer=cmat_rank[,match(drug_details_cancer$broad_id_trimmed, colnames(cmat))]
