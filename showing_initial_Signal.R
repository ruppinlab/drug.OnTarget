#Step 1 Results:: Showing that our analysis has a signal
# Total genes
total=length(df2plot$corrRank_min <1833)
# Drugs with known target as the predicted target
sum(df2plot$corrRank_min ==1)/total
# Drugs with known target as the top 10% of the ranked hits 
sum(df2plot$corrRank_min<1833)/total
# Drugs with known target as the top 20% of the ranked hits 
sum(df2plot$corrRank_min<1833*2)/total
# Drugs with known target corr >0.25
sum(figure2b_crispr$Score>0.25)

