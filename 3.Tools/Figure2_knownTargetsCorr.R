library("grid")
require('ggrepel')
source('/Users/sinhas8/myCustom_functions.R')
Targets_list = lapply(as.character(onTarget$Annotated_Target),
                      function(x) gsub(' ','',strsplit(x, '\\,')[[1]]) )
###########################################################################
# Step 1B: Known_Target_Correlation - CRISPR
###########################################################################
df2plot=onTarget$Annotated_Target_corrMeasures
head(onTarget$Annotated_Target_corrMeasures)
df2plot=df2plot[!is.infinite(df2plot$corr_max),]
df2plot=df2plot[order(df2plot$corr_max, decreasing = T),]

df2plot=df2plot[df2plot$drugCategory=='targeted cancer',]
Up_Threshold=0.25
###########################################################################
# Step 1C: Known_Target_Correlation - shRNA
###########################################################################
# For shRNA
df2plot_shRNA=onTarget$Annotated_Target_corrMeasures_shRNA
df2plot_shRNA=df2plot_shRNA[!is.infinite(df2plot$corr_max),]
df2plot_shRNA=df2plot_shRNA[order(df2plot_shRNA$corr_max, decreasing = T),]
df2plot_shRNA=df2plot_shRNA[df2plot_shRNA$drugCategory=='targeted cancer',]

Up_Threshold=0.25
###########################################################################
# Individual Figures v1
###########################################################################
# plot1<-ggplot(df2plot_shRNA, aes(y=corr_max, x=corrRank_min, shape= corr_max> Up_Threshold, color=drugCategory))+
#   geom_point()+
#   theme_bw(base_size = 15)+
#   labs(y='Known Target Essenttiality\n vs Drug Response Corr-Rho', x='Global Target gene Rank')+
#   geom_hline(yintercept = Up_Threshold)
# 
# tiff('/Users/sinhas8/Project_OffTarget/4.Results/Figure2A_shRNA.tiff', width = 900)
# plot1+  geom_label_repel(data = subset(df2plot_shRNA, corr_max>0.3),
#                          aes(label = Best_among_Annotated_Target), size=3.5, nudge_y = 0.2, nudge_x = 4000)
# dev.off()


# tiff('/Users/sinhas8/Project_OffTarget/4.Results/Figure2B.tiff')
# ggplot(df2plot, aes(y=Score, x=corrRank_min, color= corrRank_min == 1))+
#   geom_point()+
#   theme_bw(base_size = 15)+
#   labs(y='Top identified Target Correlation Rho', x='Known Target Rank')+
#   geom_hline(yintercept = Up_Threshold)+
#   theme(legend.position = "none")
# dev.off()
# 
# 

###########################################################################
# shRNA vs CRISPR comparison - Combined Figure 2A
###########################################################################
# For shRNA
scatterPlot<-ggplot(df2plot_shRNA, aes(y=corr_max, x=corrRank_min, 
                                       # shape= corr_max> Up_Threshold, 
                                       color=drugCategory))+
  geom_point()+
  theme_bw(base_size = 25)+
  labs(y='Known Target Essenttiality\n vs Drug Response Corr-Rho', x='Global Target gene Rank')+
  geom_hline(yintercept = Up_Threshold)+
  geom_label_repel(data = subset(df2plot_shRNA, corr_max>0.3),
                   aes(label = Best_among_Annotated_Target), size=4, nudge_y = 0.2, nudge_x = 5000)+
  theme(legend.position = "left")+
  ggtitle('shRNA-based')
ydensity_shRNA <- ggplot(df2plot_shRNA, aes(corr_max, fill=drugCategory)) + 
  geom_density(alpha=.5) +
  theme_bw(base_size = 25)+
  theme(legend.position = "none")+
  coord_flip()
  
A<-grid.arrange(scatterPlot, ydensity_shRNA, 
                ncol=2, nrow=1, widths=c(4, 1.4))
scatterPlot<-ggplot(df2plot, aes(y=corr_max, x=corrRank_min,
                                 # shape= corr_max> Up_Threshold,
                                 color=drugCategory))+
  geom_point()+
  theme_bw(base_size = 25)+
  labs(y='Known Target Essenttiality\n vs Drug Response Corr-Rho', x='Global Target gene Rank')+
  geom_hline(yintercept = Up_Threshold)+
  geom_label_repel(data = subset(df2plot, corr_max>0.3),
                   aes(label = Best_among_Annotated_Target), size=4, nudge_y = 0.2, nudge_x = 5000)+
  theme(legend.position = "left")+
  ggtitle('crispr-based')
ydensity <- ggplot(df2plot, aes(corr_max, fill=drugCategory)) + 
  geom_density(alpha=.5) + 
  theme_bw(base_size = 25)+
  theme(legend.position = "none")+
  coord_flip()
B<-grid.arrange(scatterPlot, ydensity, 
                ncol=2, nrow=1, widths=c(4, 1.4))
tiff('/Users/sinhas8/Project_OffTarget/4.Results/together_Targeted.tiff', width = 1200, height = 1200)
grid.arrange(A, B,ncol=1)
dev.off()

###########################################################################
# shRNA vs CRISPR comparison - Combined Figure 2B
###########################################################################
figure2b_shRNA=cbind(onTarget$Top_predicted_Target_shRNA, onTarget$Annotated_Target_corrMeasures_shRNA)
figure2b_shRNA=figure2b_shRNA[!is.infinite(figure2b_shRNA$corr_max),]
figure2b_shRNA=figure2b_shRNA[order(figure2b_shRNA$Score, decreasing = T),]
figure2b_shRNA$Score=as.numeric(as.character(figure2b_shRNA$Score))
figure2b_shRNA=figure2b_shRNA[figure2b_shRNA$drugCategory=='targeted cancer',]
dim(figure2b_shRNA)
# tiff('/Users/sinhas8/Project_OffTarget/4.Results/Figure2B.tiff')
A<-ggplot(figure2b_shRNA, aes(y=Score, x=corrRank_min, color= corrRank_min == 1))+
  geom_point()+
  theme_bw(base_size = 15)+
  labs(y='Top identified Target\n Correlation Rho', x='Known Target Rank')+
  geom_hline(yintercept = Up_Threshold)+
  theme(legend.position = "none")+ggtitle('shRNA')
# dev.off()

# For CRISPR

figure2b_crispr=cbind(onTarget$Top_predicted_Target, onTarget$Annotated_Target_corrMeasures)
figure2b_crispr=figure2b_crispr[!is.infinite(figure2b_crispr$corr_max),]
figure2b_crispr=figure2b_crispr[order(figure2b_crispr$Score, decreasing = T),]
figure2b_crispr$Score=as.numeric(as.character(figure2b_crispr$Score))
figure2b_crispr=figure2b_crispr[figure2b_crispr$drugCategory=='targeted cancer',]
# tiff('/Users/sinhas8/Project_OffTarget/4.Results/Figure2B.tiff')
B<-ggplot(figure2b_crispr, aes(y=Score, x=corrRank_min, color= corrRank_min == 1))+
  geom_point()+
  theme_bw(base_size = 15)+
  labs(y='Top identified Target\n Correlation Rho', x='Known Target Rank')+
  geom_hline(yintercept = Up_Threshold)+
  theme(legend.position = "none")
# dev.off()

tiff('/Users/sinhas8/Project_OffTarget/4.Results/2B_together_targeted.tiff')
grid.arrange(A, B, ncol=1)
dev.off()


