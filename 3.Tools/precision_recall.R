# Maximizing a Precision recall curve given a gold-standard
###########################################################################
# Add Cripsr & shRNA corr Score for known Targets
###########################################################################
require(ggpubr)
source('drug.OnTarget/3.Tools/drug_type_effect.R')
df2plot$crispr_corr=sapply(1:nrow(df2plot), function(x)
  err_handle(onTarget$corrMat[as.character(df2plot$Drugbank_Gene)[x], as.character(df2plot$drugName)[x]]) ) 
df2plot$shRNA_corr=sapply(1:nrow(df2plot), function(x) 
  err_handle(onTarget$corrMat_shRNA[as.character(df2plot$Drugbank_Gene)[x], as.character(df2plot$drugName)[x]]) ) 
df2plot$crispr_corr_rank=sapply(1:nrow(df2plot), function(x) 
  err_handle(onTarget$corrMat_rank[as.character(df2plot$Drugbank_Gene)[x], as.character(df2plot$drugName)[x]]) ) 
df2plot$shRNA_corr_rank=sapply(1:nrow(df2plot), function(x) 
  err_handle(onTarget$corrMat_shRNA_rank[as.character(df2plot$Drugbank_Gene)[x], as.character(df2plot$drugName)[x]]) ) 

###########################################################################
# Subset showing signal from DrugBank
###########################################################################
cond1 = df2plot$Action %in% c('inhibitor', 'inhibitor, competitive',
                              'antagonist', 'negative modulator')
cond2 = df2plot$drugCategory %in% 'targeted cancer'
cond3=df2plot$Citations_Count>3

sum(cond1 & cond2 & cond3)
drugs_forPR_try1=df2plot[cond1 & cond2 & cond3,]

Positive_Set = drugs_forPR_try1[c(1, 2, 4)]
Positive_Set_DrugGene_Pairs=Positive_Set[,1:2]
Shuffled_Negative_Set=data.frame(Drugbank_Gene=Positive_Set_DrugGene_Pairs[sample(1:nrow(Positive_Set_DrugGene_Pairs), 100, replace = T), 1],
                                 Chemical.Name=Positive_Set_DrugGene_Pairs[sample(1:nrow(Positive_Set_DrugGene_Pairs), 100, replace = T), 2])
Shuffled_Negative_Set=unique(Shuffled_Negative_Set)
Shuffled_Negative_Set_unique=Shuffled_Negative_Set[!apply(Shuffled_Negative_Set, 1, function(x) paste(as.character(x), collapse = '_')) %in%
                                                     apply(Positive_Set_DrugGene_Pairs, 1, function(x) paste(as.character(x), collapse = '_')),]
Negative_Set_DrugGene_Pairs=Shuffled_Negative_Set_unique[
  sample(1:nrow(Shuffled_Negative_Set_unique), nrow(Positive_Set_DrugGene_Pairs)), ]
Positive_Set_DrugGene_Pairs$Label=1
Negative_Set_DrugGene_Pairs$Label=0
Complete_Set=rbind(Positive_Set_DrugGene_Pairs, Negative_Set_DrugGene_Pairs)
Complete_Set$drug_BroadID=df2plot$drugName[match(Complete_Set$Chemical.Name, df2plot$Chemical.Name)]
head(Complete_Set)
dim(onTarget$corrMat_rank)

Complete_Set$crispr_corr=sapply(1:nrow(Complete_Set), function(x) 
  err_handle(onTarget$corrMat_rank[as.character(Complete_Set$Drugbank_Gene)[x], as.character(Complete_Set$drug_BroadID)[x]]))

Complete_Set$shRNA_corr=sapply(1:nrow(Complete_Set), function(x) 
  err_handle(onTarget$corrMat_shRNA_rank[as.character(Complete_Set$Drugbank_Gene)[x], as.character(Complete_Set$drug_BroadID)[x]]))

###########################################################################
# Figure 3A - Gold-set scores
###########################################################################
crispr_score=data.frame(Complete_Set[,-6], Screen='CRISPR')
colnames(crispr_score)[5]='Corr_Rho'
shRNA_score=data.frame(Complete_Set[,-5], Screen='shRNA')
colnames(shRNA_score)[5]='Corr_Rho'
df2plot_forPlot=rbind(crispr_score, shRNA_score)
df2plot_forPlot$Label= factor(df2plot_forPlot$Label, labels = c('Negative', 'Positive'))

tiff('/Users/sinhas8/Project_OffTarget/4.Results/crisprVSshRNA_prediction_Figure3A.tiff')
ggplot(df2plot_forPlot, aes(y=Corr_Rho, fill=factor(Label), x=factor(Label)))+
  geom_boxplot()+
  facet_wrap(~Screen)+
  stat_compare_means(method='wilcox')+
  theme_bw(base_size = 25)+
  labs(y='Correlation Rho', x='Ground Truth', fill='Ground Truth')+
  theme(legend.position = 'top')
dev.off()

require(pROC)
A<-ggroc(roc(crispr_score$Label, crispr_score$Corr_Rho))+
  theme_bw(base_size = 20)+
  ggtitle('CRISPR-Based; AUC=0.73')
B<-ggroc(roc(shRNA_score$Label, shRNA_score$Corr_Rho))+
  theme_bw(base_size = 20)+
  ggtitle('shRNA-Based; AUC=0.68')
tiff('/Users/sinhas8/Project_OffTarget/4.Results/crisprVSshRNA_AUC_Figure3B.tiff')
grid.arrange(A, B)
dev.off()
