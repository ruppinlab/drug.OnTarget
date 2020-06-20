# Finding a few drugs of interest
figure2b_crispr=figure2b_crispr[match(figure2b_shRNA$drugName, figure2b_crispr$drugName),]
imprvement_Index_crispr=figure2b_crispr$Score - figure2b_crispr$corr_max
figure2b_crispr_order=figure2b_crispr[order(imprvement_Index_crispr, decreasing = T),]
imprvement_Index_shRNA=figure2b_shRNA$Score - figure2b_shRNA$corr_max
figure2b_shRNA_order=figure2b_shRNA[order(imprvement_Index_shRNA, decreasing = T),]

predTargets_bothScreens=data.frame(Drugname=onTarget$drugCategory$name[match(figure2b_crispr_order$drugName,
                                                                             onTarget$drugCategory$broad_id_trimmed)],
                                   PredTarget_crispr=figure2b_crispr_order$PredTarget,
                                   PredTarget_shRNA=figure2b_shRNA_order$PredTarget)
onTarget$drugCategory[onTarget$drugCategory$broad_id_trimmed=='K70301465',]

