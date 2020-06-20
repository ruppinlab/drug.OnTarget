# Visualize a large matrix
library(ComplexHeatmap)
tiff('/Users/sinhas8/Project_OffTarget/test1_corrmat.tiff')
Heatmap(onTarget$corrMat, name = "corr")
dev.off()