library(exprso)
suppressWarnings(RNGversion("3.5.0"))

###########################################################
### Check plMonteCarlo

load(file.path("data.RData"))

array@annot$defineCase[1:2] <- "Case"
array@annot$defineCase[25:30] <- "Control"

# Perform bootstrapping with plMonteCarlo
set.seed(12345)
ss <- ctrlSplitSet(func = "splitStratify", percent.include = 50, colBy = NULL)
fs <- ctrlFeatureSelect(func = "fsStats", top = 0, how = "t.test")
gs <- ctrlGridSearch(func = "plGrid", how = "buildLDA", top = 2, method = "mle")
boot <- plMonteCarlo(array, B = 20, ctrlSS = ss, ctrlFS = fs, ctrlGS = gs)

# Repeat bootstrapping manually
set.seed(12345)
aucs <- vector("numeric", 20)
for(b in 1:20){

  arrays.b <- splitStratify(array, percent.include = 50, colBy = NULL)
  array.b.train <- fsStats(arrays.b[[1]], top = 0, how = "t.test")
  array.b.test <- arrays.b[[2]]
  pl.b <- plGrid(array.b.train, array.b.test, how = "buildLDA", top = 2, method = "mle")
  aucs[b] <- pl.b$valid.auc
}

test_that("plMonteCarlo is grossly intact", {

  expect_equal(
    calcMonteCarlo(boot, colBy = "valid.auc"),
    mean(aucs)
  )
})

# Check calcMonteCarlo with contrived example
set.seed(12345)
ss <- ctrlSplitSet(func = "splitStratify", percent.include = 50, colBy = NULL)
fs <- ctrlFeatureSelect(func = "fsStats", top = 0, how = "t.test")
gs <- ctrlGridSearch(func = "plGrid", how = "buildLDA", top = c(4, 3, 2), method = "mle", fold = 0)
boot <- plMonteCarlo(array, B = 1, ctrlSS = ss, ctrlFS = fs, ctrlGS = gs)

test_that("plMonteCarlo returns correctly sized @machs", {

  expect_equal(
    length(boot@machs[[1]]@preFilter[[2]]),
    4
  )

  expect_equal(
    length(boot@machs[[2]]@preFilter[[2]]),
    3
  )

  expect_equal(
    length(boot@machs[[3]]@preFilter[[2]]),
    2
  )
})

test_that("calcMonteCarlo picks best CV", {

  expect_equal(
    round(calcMonteCarlo(boot, colBy = "valid.auc"), 7),
    0.7619048
  )
})

###########################################################
### Check plNested

array@annot$defineCase[1:2] <- "Case"
array@annot$defineCase[25:30] <- "Control"

# Perform cross-validation with plNested
set.seed(12345)
fs <- ctrlFeatureSelect(func = "fsStats", top = 0)
gs <- ctrlGridSearch(func = "plGrid", how = "buildLDA", top = 0, fold = NULL, method = "mle")
nest <- plNested(array, fold = 10, ctrlFS = fs, ctrlGS = gs)

# Perform cross-validation with plCV
set.seed(12345)
cv <- plCV(array, top = 0, fold = 10, how = "buildLDA", method = "mle")

test_that("plNested without fs matches plCV", {

  expect_equal(
    mean(nest$valid.acc),
    cv
  )
})

# Check calcNested with contrived example
set.seed(12345)
fs <- ctrlFeatureSelect(func = "fsStats", top = 0)
gs <- ctrlGridSearch(func = "plGrid", how = "buildSVM", top = 2, fold = 10,
                     kernel = "linear", cost = 10^(c(-10, 1)))
nest <- plNested(array, fold = 10, ctrlFS = fs, ctrlGS = gs)

test_that("calcMonteCarlo picks best CV", {

  expect_equal(
    calcNested(nest, colBy = "train.plCV"),
    mean(nest$train.plCV[seq(2, 20, 2)])
  )
})
