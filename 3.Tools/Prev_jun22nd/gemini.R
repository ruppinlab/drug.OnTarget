# Gemini vulnerabilities
sum(unlist(prob$scna[prob$genes=='TP53',])< -0.6 & unlist(prob$scna[prob$genes=='TP53',])> -1.3)/
  length(unlist(prob$scna[prob$genes=='TP53',])< -0.6 & unlist(prob$scna[prob$genes=='TP53',])> -1.3)

tcga_arm=readxl::read_xlsx('/Users/sinhas8/Downloads/Arm_Level_TCGA.xlsx')
colnames(tcga_arm)
arm_fre=apply(tcga_arm[,14:52], 2, table)
arm_fre=arm_fre[,order(arm_fre[1,])]

# 17p
gemini=readxl::read_xlsx('/Users/sinhas8/Downloads/gemini_vulnerabilities.xlsx')

gemini_17=gemini[gemini$Chromosome==17,]
write.csv(gemini_17, '/Users/sinhas8/Downloads/gemini_vulerabilites_on_Chr17.csv')
length(gemini_17$`GEMINI gene`)
table(gemini$Chromosome)[order(table(gemini$Chromosome))]
