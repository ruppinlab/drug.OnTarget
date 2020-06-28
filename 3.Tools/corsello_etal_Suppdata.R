# Integrate info from Corsello et al
bimodal_coef=readxl::read_xlsx('2.Data/43018_2019_18_MOESM2_ESM (4).xlsx', sheet = 3, skip = 2)
biomarker=readxl::read_xlsx('2.Data/43018_2019_18_MOESM2_ESM (4).xlsx', sheet = 4, skip = 2)
biomarker_curated=readxl::read_xlsx('2.Data/43018_2019_18_MOESM2_ESM (4).xlsx', sheet = 5, skip = 2)
crisprvsDrug=readxl::read_xlsx('2.Data/43018_2019_18_MOESM2_ESM (4).xlsx', sheet = 6, skip = 2)

# TEst in crisprvsDrug
crisprvsDrug=crisprvsDrug[order(abs(crisprvsDrug$correlation), decreasing = T),]
head(crisprvsDrug, 20)
drugsofInterest