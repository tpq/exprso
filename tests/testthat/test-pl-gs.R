library(exprso)
suppressWarnings(RNGversion("3.5.0"))

###########################################################
### Check plGrid

load(file.path("data.RData"))

arrays <- splitStratify(array, percent.include = 50, colBy = NULL)
array.train <- arrays[[1]]
array.test <- arrays[[2]]

pl <- plGrid(array.train,
             array.test,
             top = 0,
             how = "buildLDA",
             fold = NULL,
             method = c("mle",
                        "mve")
)

test_that("individual @machs match expectations", {

  mach <- buildLDA(array.train, top = 0, method = "mle")
  expect_equal(
    predict(pl@machs[[1]], array.test),
    predict(mach, array.test)
  )

  mach <- buildLDA(array.train, top = 0, method = "mve")
  expect_equal(
    predict(pl@machs[[2]], array.test),
    predict(mach, array.test)
  )
})

test_that("@summary match expectations", {

  mach <- buildLDA(array.train, top = 0, method = "mle")
  expect_equal(
    matrix(pl@summary[1, c("valid.acc", "valid.sens", "valid.spec", "valid.auc")]),
    matrix(calcStats(predict(mach, array.test)))
  )

  mach <- buildLDA(array.train, top = 0, method = "mve")
  expect_equal(
    matrix(pl@summary[2, c("valid.acc", "valid.sens", "valid.spec", "valid.auc")]),
    matrix(calcStats(predict(mach, array.test)))
  )
})

###########################################################
### Check plCV

set.seed(12345)

arrays <- splitStratify(array, percent.include = 100, colBy = NULL)
array.train <- arrays[[1]]
array.test <- arrays[[2]]

acc <- plCV(array.train, top = 0, how = "buildLDA", fold = 10, method = "mle")

array.train@annot$defineCase[c(1:2)] <- "Case"
array.train@annot$defineCase[c(19:20)] <- "Control"

acc.off <- plCV(array.train, top = 0, how = "buildLDA", fold = 0, method = "mle")
pl <- plGrid(array.train,
             array.test,
             top = 0,
             how = "buildLDA",
             fold = 0,
             method = "mle"
)

test_that("plCV is grossly intact", {

  expect_equal(
    acc,
    1
  )

  expect_equal(
    acc.off,
    .8
  )

  expect_equal(
    pl@summary$train.plCV,
    .8
  )
})
