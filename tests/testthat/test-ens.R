library(exprso)
suppressWarnings(RNGversion("3.5.0"))

###########################################################
### Check plMonteCarlo

load(file.path("data.RData"))

array@annot$defineCase[1:2] <- "Case"
array@annot$defineCase[29:30] <- "Control"

set.seed(12345)
arrays <- splitStratify(array, percent.include = 50, colBy = NULL)

# Build ensemble with plGrid
set.seed(12345)
pl <- plGrid(arrays[[1]], top = 0, how = "buildSVM", fold = NULL,
             kernel = "linear", cost = 10^(c(-3, 1, 3)))
ens1 <- buildEnsemble(pl)

# Build ensemble manually
set.seed(12345)
mach1 <- buildSVM(arrays[[1]], top = 0, kernel = "linear", cost = 10^-3)
mach2 <- buildSVM(arrays[[1]], top = 0, kernel = "linear", cost = 10^1)
mach3 <- buildSVM(arrays[[1]], top = 0, kernel = "linear", cost = 10^3)
ens2 <- buildEnsemble(mach1, mach2, mach3, "notMachine")

test_that("buildEnsemble and conjoin work the same", {

  expect_equal(
    conjoin(mach1, mach2, mach3, "NOTmachine"),
    ens2
  )
})

test_that("buildEnsemble methods for pl and model match", {

  expect_equal(
    predict(ens1, arrays[[2]]),
    predict(ens2, arrays[[2]])
  )
})

ens3 <- buildEnsemble(pl, colBy = "train.acc", how = .75)
ens4 <- buildEnsemble(mach3, mach2)

test_that("buildEnsemble calls pipeFilter correctly", {

  expect_equal(
    predict(ens3, arrays[[2]]),
    predict(ens4, arrays[[2]])
  )
})
