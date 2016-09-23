library(exprso)
context("build")

###########################################################
### Build and predict ExprsMachine objects

load(file.path("data.RData"))

set.seed(12345)

arrays <- splitStratify(array, percent.include = 50, colBy = NULL)
array.train <- fsStats(arrays[[1]], top = 0)
array.test <- arrays[[2]]

mach <- buildSVM(array.train, top = 0, kernel = "linear", cost = 1)
pred.train <- predict(mach, array.train)
pred.test <- predict(mach, array.test)

test_that("ExprsBinary models predict only on ExprsBinary datasets", {

  expect_error(
    predict(mach, arrayMulti)
  )
})

test_that("ExprsPredict slots make sense", {

  expect_equal(
    as.character(pred.train@pred),
    array.train$defineCase
  )

  expect_equal(
    as.character(pred.test@pred),
    array.test$defineCase
  )

  expect_equal(
    ifelse(as.vector(pred.test@decision.values) > 0, "Case", "Control"),
    array.test$defineCase
  )

  expect_equal(
    ifelse(as.vector(pred.test@probability[,"Case"]) > .5, "Case", "Control"),
    array.test$defineCase
  )
})

test_that("built models can predict correct classes", {

  set.seed(12345)
  mach <- buildNB(array.train, top = 2)
  expect_equal(
    as.character(predict(mach, array.test)@pred),
    array.test$defineCase
  )

  set.seed(12345)
  mach <- buildLDA(array.train, top = 2)
  expect_equal(
    as.character(predict(mach, array.test)@pred),
    array.test$defineCase
  )

  set.seed(12345)
  mach <- buildSVM(array.train, top = 2, kernel = "linear", cost = 1)
  expect_equal(
    as.character(predict(mach, array.test)@pred),
    array.test$defineCase
  )

  set.seed(12345)
  mach <- buildANN(array.train, top = 2, size = 3, decay = 1)
  expect_equal(
    as.character(predict(mach, array.test)@pred),
    array.test$defineCase
  )
})

array.test@annot$defineCase[1] <- "Case"
array.test@annot$defineCase[10] <- "Control"

test_that("built models can detect wrong classes", {

  set.seed(12345)
  mach <- buildNB(array.train, top = 2)
  expect_equal(
    calcStats(predict(mach, array.test), aucSkip = TRUE)$acc,
    .9
  )

  set.seed(12345)
  mach <- buildLDA(array.train, top = 2)
  expect_equal(
    calcStats(predict(mach, array.test), aucSkip = TRUE)$acc,
    .9
  )

  set.seed(12345)
  mach <- buildSVM(array.train, top = 2, kernel = "linear", cost = 1)
  expect_equal(
    calcStats(predict(mach, array.test), aucSkip = TRUE)$acc,
    .9
  )

  set.seed(12345)
  mach <- buildANN(array.train, top = 2, size = 3, decay = 1)
  expect_equal(
    calcStats(predict(mach, array.test), aucSkip = TRUE)$acc,
    .9
  )
})

###########################################################
### Build and predict ExprsModule objects

set.seed(12345)

arraysMulti <- splitStratify(arrayMulti, percent.include = 80, colBy = NULL)
arrayMulti.train <- fsSample(arraysMulti[[1]], top = 0)
arrayMulti.test <- arraysMulti[[2]]

set.seed(12345)

mach <- buildSVM(arrayMulti.train, top = 0, kernel = "linear", cost = 1)
pred.train <- predict(mach, arrayMulti.train)
pred.test <- predict(mach, arrayMulti.test)

test_that("ExprsMulti models predict only on ExprsMulti datasets", {

  expect_error(
    predict(mach, array.train)
  )
})

mach.multi <- doMulti(arrayMulti.train, top = 0, method = "buildSVM")

test_that("doMulti performs 1 vs. all build", {

  i <- 1
  array.i <- arrayMulti.train
  array.i@annot$defineCase <- ifelse(as.numeric(array.i@annot$defineCase) == i, "Case", "Control")
  class(array.i) <- "ExprsBinary"
  expect_equal(
    calcStats(predict(mach.multi[[i]], array.i))$acc,
    1
  )

  i <- 2
  array.i <- arrayMulti.train
  array.i@annot$defineCase <- ifelse(as.numeric(array.i@annot$defineCase) == i, "Case", "Control")
  class(array.i) <- "ExprsBinary"
  expect_equal(
    calcStats(predict(mach.multi[[i]], array.i))$acc,
    1
  )

  i <- 3
  array.i <- arrayMulti.train
  array.i@annot$defineCase <- ifelse(as.numeric(array.i@annot$defineCase) == i, "Case", "Control")
  class(array.i) <- "ExprsBinary"
  expect_equal(
    calcStats(predict(mach.multi[[i]], array.i))$acc,
    1
  )
})

test_that("ExprsMulti build and predict is grossly intact", {

  expect_equal(
    calcStats(predict(mach, arrayMulti.train))$acc,
    1
  )

  expect_equal(
    calcStats(predict(mach, arrayMulti.test))$acc,
    1
  )
})
