library(exprso)
context("plGridMulti")

###########################################################
### Test plGridMulti

load(file.path("data.RData"))

set.seed(1)
splitSets <- splitStratify(arrayMulti, percent.include = 67)
fs.inner <- ctrlFeatureSelect(func = "fsStats", top = 0, how = "t.test")
pl <- plGridMulti(trainingSet(splitSets), testSet(splitSets), top = c(2, 3, 0), ctrlFS = fs.inner,
                  how = "buildSVM", kernel = c("linear", "radial"), gamma = c(.1, .2))

test_that("each ExprsMachine from plGridMulti has its own fs history", {

  expect_error(
    plGridMulti(array)
  )

  expect_equal(
    as.character(class(pl@machs[[1]])),
    "ExprsModule"
  )

  expect_equal(
    as.character(class(pl@machs[[1]]@mach[[1]])),
    "ExprsMachine"
  )

  expect_equal(
    pl@machs[[1]]@mach[[1]]@preFilter[[2]],
    c("feat4", "feat1")
  )

  expect_equal(
    pl@machs[[1]]@mach[[2]]@preFilter[[2]],
    c("feat4", "feat2")
  )

  expect_equal(
    mean(pl@summary$valid.acc),
    0.9356725
  )
})

ss <- ctrlSplitSet(func = "splitStratify", percent.include = 67)
fs <- ctrlFeatureSelect(func = "fsNULL", top = 0)
gs <- ctrlGridSearch(func = "plGridMulti", how = "buildSVM", top = c(2, 3, 0), ctrlFS = fs.inner,
                     kernel = c("linear", "radial"), gamma = c(.1, .2))

set.seed(1)
splitSets <- splitStratify(arrayMulti, percent.include = 67)
pl <- plGridMulti(trainingSet(splitSets), testSet(splitSets), top = c(2, 3, 0), ctrlFS = fs.inner,
                  how = "buildSVM", kernel = c("linear", "radial"), gamma = c(.1, .2))

set.seed(1)
pmc <- plMonteCarlo(arrayMulti, B = 1, ctrlSS = ss, ctrlFS = fs, ctrlGS = gs)

test_that("plGridMulti works within plMonteCarlo", {

  expect_equal(
    pl@summary[, c("kernel", "gamma", "cost", "train.acc", "valid.acc")],
    pmc@summary[, c("kernel", "gamma", "cost", "train.acc", "valid.acc")]
  )
})

set.seed(1)
nest <- plNested(arrayMulti, ctrlFS = fs, ctrlGS = gs, fold = 10)

test_that("plGridMulti works within plNested", {

  expect_equal(
    round(mean(nest@summary$valid.acc), 4),
    round(0.9333333, 4)
  )
})
