# Integrating shRNA and CRISPR genetic screen to create a single score (max from both)
setwd('/Users/sinhas8/Project_OffTarget/')
source('drug.OnTarget/3.Tools/Root_Step0.R')
source('drug.OnTarget/3.Tools/myCustom_functions.R')
dim(onTarget$corrMat); dim(onTarget$corrMat_shRNA)
# Unequal number of genes; For Genes common to both
# we will take the max corr out of two screen (RNAi and CRISPR)


common_genes=intersect(rownames(onTarget$corrMat), rownames(onTarget$corrMat_shRNA))
crispr_unique_genes=rownames(onTarget$corrMat)[!(rownames(onTarget$corrMat) %in% rownames(onTarget$corrMat_shRNA))]
shRNA_unique_genes=rownames(onTarget$corrMat_shRNA)[!(rownames(onTarget$corrMat_shRNA) %in% rownames(onTarget$corrMat))]
infunc_mat=rbind(onTarget$corrMat[crispr_unique_genes,],onTarget$corrMat_shRNA[shRNA_unique_genes,])
corrMat_integrated=sapply(common_genes, function(x) pmax(onTarget$corrMat[x,], onTarget$corrMat[x,]))
corrMat_integrated_t=t(corrMat_integrated)
corrMat_bothScreens=rbind(corrMat_integrated_t, infunc_mat)
corrMat_bothScreens=corrMat_bothScreens[order(rownames(corrMat_bothScreens)),]
onTarget$corrMat_bothScreens=corrMat_bothScreens
