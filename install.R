# Before using Exprso, you must install the following packages:

source("https://bioconductor.org/biocLite.R")
biocLite("affy")
biocLite("Biobase")
biocLite("limma")

install.packages(c("plyr", "lattice", "sampling", "Matching", "penalizedSVM", "pathClass", "mRMRe", "e1071", "nnet", "ROCR"))

# Check installed packages
require(plyr)
require(lattice)
require(Biobase)
require(sampling)
require(Matching)
require(penalizedSVM)
require(pathClass)
require(limma)
require(mRMRe)
require(e1071)
require(nnet)
require(ROCR)
