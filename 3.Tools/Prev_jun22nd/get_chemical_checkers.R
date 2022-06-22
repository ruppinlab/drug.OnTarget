# Playing Chemical Checkers
source("http://bioconductor.org/biocLite.R")
biocLite("rhdf5")
library(rhdf5)
h5ls("/Users/sinhas8/Downloads/signD1.h5")
mydata1 <- h5read("/Users/sinhas8/Downloads/signD1.h5", "/V")
myhead(mydata1)
mydata2 <- h5read("/Users/sinhas8/Downloads/signD1.h5", "/keys")
