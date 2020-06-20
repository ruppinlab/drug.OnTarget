GOI='F'
Confidence_Score=rowMeans(data.frame(onTarget$corrMat[GOI,],onTarget$corrMat_shRNA[GOI,]))
drugs_inOrder=data.frame(Drugname=colnames(onTarget$corrMat)[order(Confidence_Score, decreasing = T)],
                         Confidence_Score)
drugs_inOrder$commonName=onTarget$drugsCommonName[match(drugs_inOrder$Drugname, colnames(onTarget$corrMat))]
drugs_inOrder$commonName[grep('favip',drugs_inOrder$commonName)]

DOI='metformin'
DOI=grep('metformin',onTarget$drugsCommonName)
Confidence_Score=rowMeans(data.frame(onTarget$corrMat[,DOI],onTarget$corrMat_shRNA[,DOI]))
drugs_inOrder=data.frame(Drugname=colnames(onTarget$corrMat)[order(Confidence_Score, decreasing = T)],
                         Confidence_Score)
drugs_inOrder$commonName=onTarget$drugsCommonName[match(drugs_inOrder$Drugname, colnames(onTarget$corrMat))]
drugs_inOrder$commonName[grep('favip',drugs_inOrder$commonName)]
head(sort(onTarget$corrMat[,DOI], decreasing = T), 20)

