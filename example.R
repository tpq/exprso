# Source exprso from github
source("https://raw.githubusercontent.com/tpq/exprso/master/exprso.R")

# Import Golub ALL/AML data
require(golubEsets)
data(Golub_Merge)
array <- arrayEset(Golub_Merge, colBy = "ALL.AML", include = list("ALL", "AML"))
array <- modFilter(array, 20, 16000, 500, 5) # pre-filter Golub ala Deb 2003
array <- modTransform(array) # lg transform
array <- modNormalize(array, c(1, 2)) # normalize gene and subject vectors

# Several options for splitting data
arrays <- splitStratify(array, percent.include = 67, colBy = NULL)
splitSample(array, percent.include = 67, replace = FALSE)

# Option to balance the test set
array.test <- splitStratify(arrays[[2]], percent.include = 100, colBy = NULL)[[1]]

# Perform t-test feature selection
array.train <- fsStats(arrays[[1]], probes = 0, how = "t.test")

# PCA top 50 probes by t-test
array.train <- fsPrcomp(array.train, probes = 50)

# Visualization 2D (top 2 components)
layout(matrix(c(1), 1, 1, byrow = TRUE))
colors <- ifelse(array.train@annot$defineCase %in% "Case", "#00FFFFFF", "#FF0000FF")
plot(array.train@exprs[1, ], array.train@exprs[2, ], col = colors, pch = 20)

# Visualization 3D (top 3 components)
plot(array.train, i = 1, j = 2, k = 3)

# Call a 'build' function to build single classifier
mach <- buildSVM(array.train, probes = 5, kernel = "linear", cost = 1)
pred <- predict(mach, array.train) # predict model on training set
pred <- predict(mach, array.test) # predict model on test set
calcStats(pred, array.test) # calculate acc, sens, spec, auc

# You can also search through several parameters in order to build multiple machines
# Set fold = 10: run inner-loop of 10-fold cv to inform parameter selection
pl <- plGrid(array.train,
             array.test, # optional
             probes = c(3, 5), # classifier size n
             how = "buildSVM", # build function
             kernel = "linear", # kernels
             cost = 10^(-3:3), # costs
             fold = 10)

# Get probes used to build the i-th classifier
getProbeSet(pl, 4) # get the probes used to build the ith machine

# Two ways to build ensembles: (I) Manually, or (II) Systematically
ens <- buildEnsemble(pl@machs[[3]], pl@machs[[4]], pl@machs[[5]]) # manually
ens <- buildEnsemble(pl, top.N = 9, colBy = "train.plCV") # systematically

# Can also build ensembles based on the product of multiple columns
ens <- buildEnsemble(pl, top.N = 9, colBy = c("train.plCV", "valid.auc"))

# We can predict class membership using the ensemble
pred <- predict(ens, array.test)
calcStats(pred, array.test)

# We can perform high-throughput classification with resampling using plMonteCarlo
# plMonteCarlo calls 3 separate functions during the bootstrapping
#   1. It calls a split_ function to split data into a training subset and validation set
#   2. It calls N fs_ functions to process the data through feature selection
#   3. It calls a pl_ function to perform classification (i.e. plGrid)
# Each of these function calls get controlled by an argument handler
ss <- ctrlSplitSet(func = "splitStratify", percent.include = 67, colBy = NULL)
fs <- list(ctrlFeatureSelect(func = "fsStats", probes = 0, how = "t.test"),
           ctrlFeatureSelect(func = "fsPrcomp", probes = 50))
gs <- ctrlGridSearch(func = "plGrid", how = "buildSVM", probes = c(2, 3, 4), kernel = c("linear", "radial"),
                     cost = 10^(-3:3), gamma = 10^(-1:1), fold = 10)
boot <- plMonteCarlo(arrays[[1]], B = 3, ctrlSS = ss, ctrlFS = fs, ctrlGS = gs, save = FALSE)

# Running nested cross-validation works similarly
fs <- ctrlFeatureSelect(func = "fsEbayes", probes = 0)
gs <- ctrlGridSearch(func = "plGrid", how = "buildANN", probes = c(10, 20, 30), size = 1:3, decay = c(0, .5, 1), fold = 0)
nest <- plNested(arrays[[1]], fold = 10, ctrlFS = fs, ctrlGS = gs, save = FALSE)

# Ensembles built from plBoot include the top.N from each boot 
ens <- buildEnsemble(boot, top.N = 3, colBy = "valid.auc")
length(ens@machs) # Total of top.N = 3 * B = 3...9 machs!
pred <- predict(ens, array.test) # As above...
calcStats(pred, array.test)

# We can also filter and summarize ExprsPipeline objects
bootSubset <- pipeFilter(boot, colBy = "valid.auc", how = 75, gate = 95)
summary(bootSubset)
