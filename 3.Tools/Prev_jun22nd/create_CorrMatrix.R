# Project Off Target
###########################################################################
# Step 0: Script to Load Basic Set of libraries, Variables and 
###########################################################################
# Libraries
require(parallel); require(tictoc); require(WGCNA); require(ggplot2); require(cowplot)
source('/Users/sinhas8/myCustom_functions.R')
# # Files
setwd('/Users/sinhas8/Project_OffTarget/2.Data/')
depmap=get(load('depmap_19Q3.RData'))
crispr=depmap$avana[,intersect(colnames(depmap$avana), colnames(depmap$drug_prism))]
prism=depmap$drug_prism[,intersect(colnames(depmap$avana), colnames(depmap$drug_prism))]
###########################################################################
# Step 0: Script to Load Basic Set of libraries, Variables and Function for Project Off Target
###########################################################################
use_cores=detectCores()
# use_cores=4
crispr_list=lapply(split(crispr, c(rep(1:use_cores, each=nrow(crispr)/use_cores),
                                   rep(use_cores, nrow(crispr)%%use_cores))),
                   matrix, ncol=ncol(crispr))
calling_CORtest <- function(crispr_mat){
  lapply(1:nrow(crispr_mat), function(y) sapply(1:nrow(prism), function(x)
    cor(crispr_mat[y,],
        prism[x,],
        use = "pairwise.complete.obs",
        quick = 1,
        weights.x = NULL,
        weights.y = NULL)
  ))
}
tic()
tt=mclapply(1:length(crispr_list), function(x) err_handle(calling_CORtest(crispr_list[[x]])) ,
            mc.cores = detectCores())
# tt=calling_CORtest(crispr_list[[1]])
toc()
saveRDS(tt, '/Users/sinhas8/Project_OffTarget/2.Data/corrMat.RDS')
###########################################################################
# Step 1B: Cluster and visualize the matrix
###########################################################################
# cmat=readRDS('/Users/sinhas8/Project_OffTarget/2.Data/corrMat.RDS')
# cmat=unlist(cmat, recursive = F)
# cmat=do.call(rbind, cmat)
# saveRDS(cmat, '/Users/sinhas8/Project_OffTarget/2.Data/corrMat_matricForm.RDS')
cmat=readRDS('/Users/sinhas8/Project_OffTarget/2.Data/corrMat_matricForm.RDS')
colnames(cmat)=rownames(prism)
rownames(cmat)=rownames(crispr)
###########################################################################
# Step 1A: Cluster and visualize the matrix
###########################################################################
plotHist<-function(testVector, xtitle){
  cmat_colMaximum=data.frame(value=testVector)
  ggplot(cmat_colMaximum, aes(x=value))+
    geom_freqpoly(bins=40)+
    geom_histogram(color="black", fill="white", bins=50)+
    labs(x=xtitle, title=xtitle)
}
cmat_colMax=apply(cmat, 2, max)
cmat_colMedian=apply(cmat, 2, median)
cmat_colMean=apply(cmat, 2, mean)
cmat_colSD=apply(cmat, 2, sd)
tiff('/Users/sinhas8/Project_OffTarget/4.Results/corMat_landscape.tiff', 
     height=600, width = 600)
plot_grid(plotHist(cmat_colMedian, xtitle = 'colMedian'),
          plotHist(cmat_colMean, xtitle = 'colMean'),
          plotHist(cmat_colMax, xtitle = 'colMax'),
          plotHist(cmat_colSD, xtitle = 'colSD')
)
dev.off()
###########################################################################
# Step 1B: Target Correlation
###########################################################################
Targets_list=lapply(as.character(depmap$drug_target_map$target),
                    function(x) gsub(' ','',strsplit(x, '\\,')[[1]]) )
Annotated_Target_corr=sapply(1:length(Targets_list), function(x)
  sapply(Targets_list[[x]], function(y)  err_handle(onTarget$corrMat[y,x]) ))
Annotated_Target_corr_mean=sapply(Annotated_Target_corr, function(x) mean(x, na.rm=T))
Annotated_Target_corr_max=sapply(Annotated_Target_corr, function(x) max(x, na.rm=T))
tiff('/Users/sinhas8/Project_OffTarget/4.Results/Target_annotation_corMat_landscape.tiff', 
     height=600, width = 600)
plot_grid(plotHist(Annotated_Target_corr_mean, xtitle='Annotated_Target_corr_mean'),
          plotHist(Annotated_Target_corr_max,xtitle='Annotated_Target_corr_max'),
          plotHist(Annotated_Target_corr_rank_mean, xtitle='Annotated_Target_corr_mean_RANK'),
          plotHist(Annotated_Target_corr_rank_max,xtitle='Annotated_Target_corr_max_RANK'))
dev.off()
###########################################################################
# Step 1C: Target Corr Rank 
###########################################################################
cmat_rank=apply(cmat, 2, function(x) rank(x) )
Annotated_Target_corr_rank=sapply(1:length(Targets_list), function(x)
  sapply(Targets_list[[x]], function(y)  err_handle(cmat_rank[y,x]) ))
Annotated_Target_corr_rank_mean=sapply(Annotated_Target_corr_rank, function(x) mean(x, na.rm=T))
Annotated_Target_corr_rank_max=sapply(Annotated_Target_corr_rank, function(x) max(x, na.rm=T))
plotHist()
# tiff('/Users/sinhas8/Project_OffTarget/4.Results/Target_annotation_corRANKMat_landscape.tiff', 
#      height=600, width = 600)
# plot_grid(
# dev.off()
###########################################################################
# Perfect annotation
###########################################################################
Annotated_Target_corr[which(Annotated_Target_corr_rank_mean== max(Annotated_Target_corr_rank_mean, na.rm = T) )]
# Semi-Perfect annotation
Annotated_Target_corr[which(Annotated_Target_corr_rank_mean>max(Annotated_Target_corr_rank_mean, na.rm = T)-1000)]
# Semi-Perfect annotation
(Annotated_Target_corr[depmap$drug_target_map$phase=='Launched'])[which(Annotated_Target_corr_rank_mean[depmap$drug_target_map$phase=='Launched']>max(Annotated_Target_corr_rank_mean, na.rm = T)-200)]
###########################################################################
# Recent Remapping
###########################################################################
# SGI-1776
cmat_rank['PIM1',grep('1776', depmap$drug_target_map$name, ignore.case = T)]
# nutlin-3a to TP53
which.max(cmat_rank[,grep('nutlin-3', depmap$drug_target_map$name, ignore.case = T)])
###########################################################################
# Only working with Oncology drugs
###########################################################################
head(depmap$drug_target_map[order(Annotated_Target_corr_rank_mean, decreasing = T),]$name, 50)
sum(grepl('inhibitor',depmap$drug_target_map$moa))
sum(grepl('agonist',depmap$drug_target_map$moa))
sum(grepl('antagonist',depmap$drug_target_map$moa))

