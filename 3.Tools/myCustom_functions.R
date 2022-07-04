###Following are custom functions I've made to avoid repetitively writing same things again - after three years of my PhD :p
myhead<-function(x){
  x[1:min(5, nrow(x)), 1:min(5, ncol(x))]
}
err_handle<-function(x){ tryCatch(x, error=function(e){NA}) }
test_topbottomQuantile<-function(var1=UCD_score_M1, var2_tobetiled=unlist(Exp_metab[1,]), which_tail='g', numofQuantiles=3){
  length(var1); length(var2_tobetiled)
  var2_tobetiled=xtile(var2_tobetiled, numofQuantiles)
  ##Only keep top/bottom quantile
  var1= var1[var2_tobetiled== min(var2_tobetiled) | var2_tobetiled== max(var2_tobetiled)]
  var2_tobetiled= var2_tobetiled[var2_tobetiled== min(var2_tobetiled) | var2_tobetiled== max(var2_tobetiled)]
  c(sig=err_handle(wilcox.test(var1 ~ factor(var2_tobetiled), alternative=which_tail)$p.value),
    eff_size=err_handle(diff(aggregate(var1, by=list(var2_tobetiled), mean)[,2]) ) )
}
hypergeometric_test_for_twolists<-function(test_list, base_list, global, lowertail=FALSE) {
  #If lowertail=FALSE - we calculate the probability for enrichment
  length(base_list)
  base_list_within_global=global[na.omit(match(base_list, global))]
  intersect_of_two_list= test_list[!is.na(match(test_list, base_list_within_global))]
  phyper(length(intersect_of_two_list)-1, #white balls in the samples
         length(base_list_within_global), #total white balls in the box
         length(global)- length(base_list_within_global), #total black balls in the box
         length(test_list), #total balls sampled
         lower.tail=lowertail) #whether you wish to calculate enrichment or depletion
}
fdrcorr<-function(test_list){p.adjust(test_list, method = 'fdr')}

# dftemp=data.frame(prob$samples, prob$types, prob$stage)
# sum(aggregate(prob.samples ~ prob.stage+prob.types, data=dftemp[dftemp$prob.types=='LUSC',], function(x) length(x))[1:3,3])/sum(aggregate(prob.samples ~ prob.stage+prob.types, data=dftemp[dftemp$prob.types=='LUSC',], function(x) length(x))[,3])
# sum(aggregate(prob.samples ~ prob.stage+prob.types, data=dftemp[dftemp$prob.types=='LUAD',], function(x) length(x))[1:4,3])/sum(aggregate(prob.samples ~ prob.stage+prob.types, data=dftemp[dftemp$prob.types=='LUAD',], function(x) length(x))[,3])
#split(aggregate(mat$X ~ mat$stage_trimmed+mat$hist, data=dftemp, function(x) length(x)))

#Subsetting a set of columns
colSubset<-function(mat, column_Names){
  mat[,na.omit(match(column_Names, colnames(mat)))]
}
#Subsetting a set of rows
rowSubset<-function(mat, row_Names){
  mat[na.omit(match(row_Names, rownames(mat))),]
}

vectorSubset<-function(vec, Names){
  vec[!is.na(match(names(vec), Names))]
}


installORload<-function(packages){
  package.check <- lapply(packages, FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  })
}
# Scale 0-1
topXPercentValue<-function(vec, X_percentile=95){
  vec=na.omit(vec)
  len=length(vec)
  vec=sort(vec)
  vec[ceiling(len*(X_percentile/100))]
}
range01 <- function(x){
  # Chossing 95% and 5% percentile as thresholds for outliers
  substitute_of_Min=topXPercentValue(vec=x, 
                                     X_percentile=5)
  substitute_of_Max=topXPercentValue(vec=x, 
                                     X_percentile=95)
  x_scaled=(x-substitute_of_Min)/(substitute_of_Max-substitute_of_Min)
  x_scaled[x_scaled<0]=0
  x_scaled[x_scaled>1]=1
  x_scaled
}

range01(1:200)


# Scale 0-1 cancer type specifically
scaling_cancerType<-function(quan1=gi, quan2=hist){
  unlist(lapply(split(quan1, quan2), function(x) range01(x)))
}


factor2numeric<-function(x){
  as.numeric(as.character(x))
}

mycorTest<-function(x, y){
  unlist(cor.test(x, y)[c(3, 4)])
}

stripall2match<-function(x){
  # Strip all non-char and non-numeric and make lower case
  # this is primarily to facilitate inconsistent naming
  tolower(gsub('[^A-z0-9]','',x) )
}


'%!in%' <- function(x,y)!('%in%'(x,y))


strsplit_customv0 <- function(infunc_list=pred_viab$cellLines_mapping$cellLine_ID,
                              infunc_split_by='_',
                              retreving_onject_id=1){
  sapply(strsplit(infunc_list, split = infunc_split_by), function(x) x[retreving_onject_id])
}

colMax <- function (colData) {
  apply(colData, MARGIN=c(2), max)
}
colMedian <- function (colData) {
  apply(colData, MARGIN=c(2), median)
}

rowMax <- function (colData) {
  apply(colData, MARGIN=c(1), max)
}
rowMin <- function (colData) {
  apply(colData, MARGIN=c(1), min)
}


srowProd <- function (colData) {
  apply(colData, 1, prod)
}


match_stripall <- function(x, y){
  x_striped=tolower(gsub('[^A-z0-9]','',x) )
  y_striped=tolower(gsub('[^A-z0-9]','',y) )
  match(x_striped, y_striped)
}
