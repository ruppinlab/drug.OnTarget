# Integrating shRNA and CRISPR genetic screen to create a single score (max from both)
setwd('/Users/sinhas8/Project_OffTarget/3.Tools/')
source('Root_Step0.R')
source('myCustom_functions.R')
dim(onTarget$corrMat); dim(onTarget$corrMat_shRNA)
# Unequal number of genes; For Genes common to both
# we will take the max corr out of two screen (RNAi and CRISPR)

max_outofTwoScreens
common_genes=intersect(rownames(onTarget$corrMat), rownames(onTarget$corrMat_shRNA))
crispr_unique_genes=rownames(onTarget$corrMat)[!(rownames(onTarget$corrMat) %in% rownames(onTarget$corrMat_shRNA))]
shRNA_unique_genes=rownames(onTarget$corrMat_shRNA)[!(rownames(onTarget$corrMat_shRNA) %in% rownames(onTarget$corrMat))]

infunc_mat=rbind(onTarget$corrMat[crispr_unique_genes,],onTarget$corrMat_shRNA[shRNA_unique_genes,])
dim(infunc_mat)