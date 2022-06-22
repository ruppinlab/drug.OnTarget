# Whether Irutinib binds to EGFR-mutant or EGFR-WT
require("DEGreport")
common_colnames=intersect(colnames(onTarget$drug_prism), colnames(onTarget$avana))
egfr_Mut=onTarget$mutations_matrix['EGFR',]==1
egfr_Mut_mapped=egfr_Mut[common_colnames]
drug_ofInterest=onTarget$drugCategory[which(onTarget$drugCategory$broad_id_trimmed=='K70301465'),]
# background Distribution
cor.test(onTarget$drug_prism[as.character(drug_ofInterest$broad_id_trimmed),common_colnames],
         onTarget$avana['EGFR',common_colnames],method = 's')
# In context of Mutation
cor.test(onTarget$drug_prism[as.character(drug_ofInterest$broad_id_trimmed),names(which(egfr_Mut_mapped))],
         onTarget$avana['EGFR',names(which(egfr_Mut_mapped))],method = 's')
# In context of WT
cor.test(onTarget$drug_prism[as.character(drug_ofInterest$broad_id_trimmed),names(which(!egfr_Mut_mapped))],
         onTarget$avana['EGFR',names(which(!egfr_Mut_mapped))],method = 's')


df2plot=data.frame(irutinib_response=onTarget$drug_prism[as.character(drug_ofInterest$broad_id_trimmed),common_colnames],
                   egfr_ess=onTarget$avana['EGFR',common_colnames],
                   mutation_status=as.numeric(egfr_Mut_mapped))

ggplot(df2plot, aes(x=irutinib_response, y=egfr_ess))+
  geom_point()+
  stat_smooth(method = 'lm')+
  theme_bw(base_size = 15)+
  facet_grid(.~mutation_status)+
  geom_cor(method='p')

# Conclusion - EGFR mutant is stronger