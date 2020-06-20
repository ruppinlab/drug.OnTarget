# Validation V2 - Gold Standard based on pharmacologically active drugs from DrugBank
require(ggpubr)
load('/Users/sinhas8/Downloads/drugs2targets.RData')
active_drugs2targets=drugs2targets[which(drugs2targets$pharmaco.active=='yes' & 
                                           drugs2targets$action.simp.strict=='inhibition' ),]
# source('/Users/sinhas8/Project_OffTarget/3.Tools/fewDrugs_for_experiements.R')
matched_gs3=active_drugs2targets[match(tolower(predTargets_bothScreens$Drugname), tolower(active_drugs2targets$drug.name)),]
Positive_Set_DrugGene_Pairs=na.omit(data.frame(Drugbank_Gene=matched_gs3$protein.gene.symbol, 
                                               Chemical.Name=matched_gs3$drug.name))

# Positive_Set_DrugGene_Pairs=Positive_Set[,1:2]

Shuffled_Negative_Set=data.frame(Drugbank_Gene=Positive_Set_DrugGene_Pairs[sample(1:nrow(Positive_Set_DrugGene_Pairs), 100, replace = T), 1],
                                 Chemical.Name=Positive_Set_DrugGene_Pairs[sample(1:nrow(Positive_Set_DrugGene_Pairs), 100, replace = T), 2])
Shuffled_Negative_Set=unique(Shuffled_Negative_Set)

Shuffled_Negative_Set_unique=Shuffled_Negative_Set[!apply(Shuffled_Negative_Set, 1, function(x) paste(as.character(x), collapse = '_')) %in%
                                                     apply(Positive_Set_DrugGene_Pairs[,1:2], 1, function(x) paste(as.character(x), collapse = '_')),]

Negative_Set_DrugGene_Pairs=Shuffled_Negative_Set_unique[
  sample(1:nrow(Shuffled_Negative_Set_unique), nrow(Positive_Set_DrugGene_Pairs)), ]

Positive_Set_DrugGene_Pairs$Label=1
Negative_Set_DrugGene_Pairs$Label=0
Complete_Set=rbind(Positive_Set_DrugGene_Pairs, Negative_Set_DrugGene_Pairs)

Complete_Set$drug_BroadID=onTarget$drugCategory$broad_id_trimmed[match(Complete_Set$Chemical.Name, tolower(onTarget$drugCategory$name))]

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

tiff('/Users/sinhas8/Project_OffTarget/4.Results/crisprVSshRNA_prediction_Figure3A_PharmaActive_DrugBank.tiff')
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
  ggtitle('CRISPR-Based; AUC=0.66')
B<-ggroc(roc(shRNA_score$Label, shRNA_score$Corr_Rho))+
  theme_bw(base_size = 20)+
  ggtitle('shRNA-Based; AUC=0.67')
tiff('/Users/sinhas8/Project_OffTarget/4.Results/crisprVSshRNA_AUC_Figure3B_PharmaActiveCompounds.tiff')
grid.arrange(A, B)
dev.off()
