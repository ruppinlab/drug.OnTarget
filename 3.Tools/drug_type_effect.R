# What is the effect of Drug Type?
###########################################################################
# shRNA vs CRISPR - Type of Drugs - Combined Figure 2C
###########################################################################
source('drug.OnTarget/3.Tools/ValidationA1_DrugBankCitations.R')
DrugBankScore_vsOurScore=data.frame(biogrid_matched_Subsetlist_withScore, 
                                    PredvsKnown_scores_matched[match(
                                      tolower(biogrid_matched_Subsetlist_withScore$Chemical.Name),
                                      tolower(PredvsKnown_scores_matched$CommonDrugName)),])

DrugBankScore_vsOurScore$drugCategory=onTarget$drugCategory$drug_category[
  match(DrugBankScore_vsOurScore$drugName, colnames(onTarget$corrMat))]

colnames(DrugBankScore_vsOurScore)[c(1)]='Drugbank_Gene'
colnames(DrugBankScore_vsOurScore)[c(4)]='Citations_Count'

df2plot=as.data.frame(DrugBankScore_vsOurScore)
df2plot$Count_KnownTarget=sapply(df2plot$KnownTarget,
                                 function(x) length(strsplit(as.character(x), '\\,')[[1]])) 

# df2plot_subset$drugName[!duplicated(df2plot_subset$drugName)]
# df2plot_subset[df2plot_subset$drugName=='K18895904',]
###########################################################################
# Try1 - Figure 2C
###########################################################################
tiff('/Users/sinhas8/Project_OffTarget/4.Results/together.tiff', width = 1200, height = 1200)
ggplot(df2plot_subset, aes(y=KnownTarget_corrMean, x=Citations_Count %/%5 +1,
                           # color = factor(xtile(Count_KnownTarget,5))
                           color = Action
                           ))+
  geom_point()+
  facet_wrap(~drugCategory, nrow=2)+
  theme_bw(base_size = 25)+
  labs(y='Known Target Essenttiality\n vs Drug Response Corr-Rho', x='Number of Citations')+
  geom_hline(yintercept = 0.25)
  # geom_label_repel(data = subset(df2plot, KnownTarget_corrMax>0.3),
  #                  aes(label = Best_among_KnownTarget_based_onCorr, size=4, nudge_y = 0.2, nudge_x = 5000)+
  # theme(legend.position = "left")+
  # ggtitle('crispr-based')
dev.off()
###########################################################################
# Try2 - Figure 2C
###########################################################################
df2plot_subset$crispr_corr=sapply(1:nrow(df2plot_subset), function(x) 
  err_handle(onTarget$corrMat[as.character(df2plot_subset$Drugbank_Gene)[x], as.character(df2plot_subset$drugName)[x]]) ) 
df2plot_subset$crispr_corr_rank=sapply(1:nrow(df2plot_subset), function(x) 
  err_handle(onTarget$corrMat_rank[as.character(df2plot_subset$Drugbank_Gene)[x], as.character(df2plot_subset$drugName)[x]]) ) 
df2plot_subset$shRNA_corr_rank=sapply(1:nrow(df2plot_subset), function(x) 
  err_handle(onTarget$corrMat_shRNA_rank[as.character(df2plot_subset$Drugbank_Gene)[x], as.character(df2plot_subset$drugName)[x]]) ) 
Up_Threshold=0.25
tiff('/Users/sinhas8/Project_OffTarget/4.Results/together_Try2.tiff', width = 1200, height = 1200)
ggplot(df2plot_subset, aes(y = crispr_corr, x= factor(Citations_Count >5),
                           # fill = factor(Count_KnownTarget>1)
                            fill = Action %in% c('inhibitor', 'inhibitor, competitive',
                                                 'antagonist', 'negative modulator')))+
  geom_boxplot()+
  facet_wrap(~drugCategory, nrow=2)+
  theme_bw(base_size = 25)+
  labs(y='Known Target Essenttiality\n vs Drug Response Corr-Rho', x='Number of Citations')+
  geom_hline(yintercept = Up_Threshold)
# geom_label_repel(data = subset(df2plot, KnownTarget_corrMax>0.3),
#                  aes(label = Best_among_KnownTarget_based_onCorr, size=4, nudge_y = 0.2, nudge_x = 5000)+
# theme(legend.position = "left")+
# ggtitle('crispr-based')
dev.off()

###########################################################################
# Try3 - Figure 2C
###########################################################################
cond1 = df2plot$Action %in% c('inhibitor', 'inhibitor, competitive',
                              'antagonist', 'negative modulator')
cond2 = df2plot$drugCategory %in% 'targeted cancer'
df2plot$Citations_Count_discrete=factor(df2plot$Citations_Count >5, labels = c('Low', 'High'))
head(df2plot)
tiff('/Users/sinhas8/Project_OffTarget/4.Results/together_Try3.tiff', width = 600, height = 600)
ggplot(df2plot[cond1 & cond2, ], aes(y = KnownTarget_corrMean, x=Citations_Count_discrete ,
                            fill = factor(Count_KnownTarget<2)) )+
  geom_boxplot()+
  geom_point(aes(fill = factor(Count_KnownTarget<2) ))+
  facet_wrap(~factor(Count_KnownTarget<2), nrow=1)+
  theme_bw(base_size = 25)+
  labs(y='Known Target Essenttiality\n vs Drug Response Corr-Rho', x='Number of Citations')+
  geom_hline(yintercept = Up_Threshold)+
  theme(legend.position = "top")
dev.off()
saveRDS(df2plot_launched, '/Users/sinhas8/Project_COVID19/DrugsList.RDS')