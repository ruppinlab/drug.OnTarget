# Identify TRUE targets of metformin
ccle$Target_mapping$name[grep('metfo',ccle$Target_mapping$name, ignore.case = T)]
length(ccle$drug[grep('metfo',ccle$Target_mapping$name, ignore.case = T),])
metformin_resp=ccle$drug[grep('metfo',ccle$Target_mapping$name, ignore.case = T),]
annotation_subset=ccle$annotation[match(colnames(ccle$drug), ccle$annotation$depMapID),]
  write.csv(data.frame(table(annotation_subset$tcga_code)), '/Users/sinhas8/Project_OffTarget/prism_cancerType_distribution.csv')
df2plot=cbind(annotation_subset, metformin_resp)
table(df2plot$tcga_code)
df2plot_hnsc=df2plot[which(df2plot$tcga_code=='HNSC'),]

ggplot(df2plot, aes(x=Supplements, y=metformin_resp))+
  geom_boxplot()+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# Deriving from corr-matrix
MetFormain_corr_withGenes=cmat[,grep('metfo',ccle$Target_mapping$name, ignore.case = T)]
MetFormain_corr_withGenes=MetFormain_corr_withGenes[order(MetFormain_corr_withGenes, decreasing = T)]
write.table(data.frame(MetFormain_corr_withGenes), 
          '/Users/sinhas8/Project_OffTarget/MetFormain_corr_withGenes.rnk', quote=F, sep='\t')
