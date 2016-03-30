# Before using Exprso, you must install the following packages:

source("https://bioconductor.org/biocLite.R")
biocLite("affy")
biocLite("AnnotationDbi")
biocLite("GEOquery")
biocLite("Biobase")
biocLite("limma")

install.packages(
  c("plyr",
    "lattice",
    "sampling",
    "Matching",
    "penalizedSVM",
    "pathClass",
    "mRMRe",
    "e1071",
    "MASS",
    "nnet",
    "randomForest",
    "ROCR")
)
