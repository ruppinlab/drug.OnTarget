# Ibrutibin_Target
# Testing the symmetric index of drug-target identification
df2plot2=data.frame(exp=onTarget$expression['BTK',intersect(colnames(onTarget$avana),
                                    colnames(onTarget$drug_prism))],
           tissue=onTarget$annotation$Site_Primary[match(intersect(colnames(onTarget$avana),
                                                              colnames(onTarget$drug_prism)),
                onTarget$annotation$depMapID)])
levels(df2plot2$tissue)= paste(names(table(df2plot2$tissue)), table(df2plot2$tissue))
tiff('4.Results/BTK_across_tissue_DOTCellLines.tiff', width=900)
ggplot(df2plot2, aes(y=exp, x=tissue))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
dev.off()
# Conclusion - BTK is only expressed in haematopeitic cell lines

# Find Irutibin target in Haep vs non-Haep cell lines
with_BTK_expression = colnames(onTarget$expression)[onTarget$expression['BTK',]>1]
COI=intersect(colnames(onTarget$avana), colnames(onTarget$drug_prism))
COI_with_BTK=intersect(COI, with_BTK_expression)
COI_without_BTK=COI[!(COI %in% with_BTK_expression)]
length(COI_without_BTK)
DOI=grep('ibrutinib',tolower(onTarget$drugCategory$name))
GOI=match(c('BTK', 'EGFR'),rownames(onTarget$avana))  

# Without BTK
in_corr_without_BTK=sapply(1:nrow(onTarget$avana), function(x) 
  unlist(cor.test(unlist(onTarget$drug_prism[DOI,COI_without_BTK]),
           unlist(onTarget$avana[x,COI_without_BTK]), method='s')[c(3,4)]))
in_corr_without_BTK=t(in_corr_without_BTK)
rownames(in_corr_without_BTK)=rownames(onTarget$avana)
# With BTK

in_corr_with_BTK=sapply(1:nrow(onTarget$avana), function(x) 
  unlist(cor.test(unlist(onTarget$drug_prism[DOI,COI_with_BTK]),
                  unlist(onTarget$avana[x,COI_with_BTK]), method='s')[c(3,4)]))
in_corr_with_BTK=t(in_corr_with_BTK)
rownames(in_corr_with_BTK)=rownames(onTarget$avana)

data.frame(in_corr_with_BTK[GOI,], in_corr_without_BTK[GOI,])


match('BTK',rownames(in_corr_with_BTK[order(in_corr_with_BTK[,2], decreasing = T),]))
match('BTK',rownames(in_corr_without_BTK[order(in_corr_without_BTK[,2], decreasing = T),]))

