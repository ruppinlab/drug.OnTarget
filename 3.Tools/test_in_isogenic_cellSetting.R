# Using isogenic cell lines
###########################################################################
# Load Dataset
###########################################################################
isoscreen=readxl::read_xlsx('/Users/sinhas8/Downloads/ScienceDirect_files_02Mar2020_19-04-52.262/1-s2.0-S2211124719306631-mmc2.xlsx',
                            skip = 1)
source('/Users/sinhas8/Project_OffTarget/3.Tools/Root_Step0.R')
source('/Users/sinhas8/myCustom_functions.R')
###########################################################################
# 
###########################################################################
match(tolower(colnames(isoscreen)), tolower(onTarget$drugsCommonName))
sapply(c('sorafenib', 'regorafenib', 'lenvatinib', 'quizartinib', 'thioguanine', 'gemcitabine', 'fluorouracil',
         'irinotecan', 'oxaliplatin', 'sch772984'), 
       function(x) grep(x, tolower(onTarget$drugsCommonName)))


