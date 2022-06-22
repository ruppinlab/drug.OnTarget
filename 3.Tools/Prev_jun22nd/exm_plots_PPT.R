# Plotting Example plots for ppt
onTarget=readRDS('/Users/sinhas8/Project_OffTarget/2.Data/onTarget_v7.RDS')
common_cellLines=intersect(colnames(onTarget$avana),
                           colnames(onTarget$drug_prism))
onTarget$drugsCommonName
df2plot=data.frame(MDM2_KO_VFC=onTarget$avana['MDM2',common_cellLines],
           idasnutlin=onTarget$drug_prism[grep('nutlin',onTarget$drugsCommonName)[2],
                                          common_cellLines],
           geneName='MDM2',
           drugName='Idasnutlin')

df2plot$cellLine_Name=''
df2plot$cellLine_Name[match(head(rownames(df2plot)[order(df2plot$idasnutlin,decreasing = F)], 10),
                            rownames(df2plot))]=head(rownames(df2plot)[order(df2plot$idasnutlin,decreasing = F)], 10)

tiff('/Users/sinhas8/Desktop/p1_mdm2.tiff')
ggplot(df2plot, aes(y=MDM2_KO_VFC,
                    x= reorder(rownames(df2plot),
                               MDM2_KO_VFC),
                    color=MDM2_KO_VFC,
                    label=cellLine_Name))+
  geom_bar(stat = "identity")+
  labs(y='Viability fold Change after knockout',x='CellLines')+
  geom_vline(xintercept = -1, linetype='dashed', color='blue')+
  theme_classic(base_size = 23)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  geom_text_repel(nudge_x = 100)
dev.off()

tiff('/Users/sinhas8/Desktop/p1_Idasnutlin.tiff')
ggplot(df2plot, aes(y=idasnutlin,
                    x= reorder(rownames(df2plot), idasnutlin),
                    color=idasnutlin,
                    label=cellLine_Name))+
  geom_bar(stat = "identity")+
  labs(y='Viability fold Change after drug treatment',x='CellLines')+
  geom_vline(xintercept = -1, linetype='dashed', color='blue')+
  theme_classic(base_size = 23)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  geom_text_repel(nudge_x = 100)
dev.off()

tiff('/Users/sinhas8/Desktop/together.tiff')
ggplot(df2plot, aes(y= idasnutlin, x= MDM2_KO_VFC))+
  geom_point()+
  stat_smooth(method='lm')+
  labs(y='VFC after drug treatment',x='VFC after CRISPR-Cas9')+
  theme_classic(base_size = 23)+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  stat_poly_eq(formula = 'y~x')
dev.off()

KO_VFC=onTarget$avana[,common_cellLines]
nutlin_resp=onTarget$drug_prism[grep('nutlin',onTarget$drugsCommonName)[2],
                                                  common_cellLines]

myhead(nutlin_resp)
all_cor_by_nutlin=apply(KO_VFC, 1, function(x) unlist(cor.test_trimmed_v0(x, nutlin_resp)))
all_cor_by_nutlin=t(all_cor_by_nutlin)
all_cor_by_nutlin=all_cor_by_nutlin[order(all_cor_by_nutlin[,2], decreasing = T),]
head(all_cor_by_nutlin[,2])
