library(exprso)
context("build")

###########################################################
### Prepare ExprsBinary and ExprsMulti objects

set.seed(1235)

df.a <- data.frame(
  "id" = 1:10,
  "class" = rep("a", 10),
  "sex" = c(rep("M", 5), rep("F", 5)),
  "feat1" = rnorm(10, mean = 10, sd = 1),
  "feat2" = rnorm(10, mean = 20, sd = 5),
  "feat3" = rnorm(10, mean = 5, sd = 1),
  "feat4" = rnorm(10, mean = 40, sd = 1)
)

df.b <- data.frame(
  "id" = 11:30,
  "class" = rep("b", 20),
  "sex" = c(rep("M", 10), rep("F", 10)),
  "feat1" = rnorm(20, mean = 20, sd = 5),
  "feat2" = rnorm(20, mean = 10, sd = 1),
  "feat3" = rnorm(20, mean = 5, sd = 1),
  "feat4" = rnorm(20, mean = 20, sd = 1)
)

df.c <- data.frame(
  "id" = 31:40,
  "class" = rep("c", 10),
  "sex" = c(rep("M", 3), rep("F", 7)),
  "feat1" = rnorm(10, mean = 15, sd = 3),
  "feat2" = rnorm(10, mean = 15, sd = 3),
  "feat3" = rnorm(10, mean = 5, sd = 1),
  "feat4" = rnorm(10, mean = 30, sd = 1)
)

df <- do.call(rbind, list(df.a, df.b, df.c))

tempFile <- tempfile()
write.table(df, file = tempFile, sep = "\t")

array <-
  arrayRead(tempFile, probes.begin = 4, colID = "id", colBy = "class",
            include = list("a", "b"))

arrayMulti <-
  arrayRead(tempFile, probes.begin = 4, colID = "id", colBy = "class",
            include = list("a", "b", "c"))

###########################################################
### Build and predict ExprsMachine objects

set.seed(12345)

arrays <- splitStratify(array, percent.include = 50, colBy = NULL)
array.train <- fsStats(arrays[[1]], probes = 0)
array.test <- arrays[[2]]

mach <- buildSVM(array.train, probes = 0, kernel = "linear", cost = 1)
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
  mach <- buildNB(array.train, probes = 2)
  expect_equal(
    as.character(predict(mach, array.test)@pred),
    array.test$defineCase
  )

  set.seed(12345)
  mach <- buildLDA(array.train, probes = 2)
  expect_equal(
    as.character(predict(mach, array.test)@pred),
    array.test$defineCase
  )

  set.seed(12345)
  mach <- buildSVM(array.train, probes = 2, kernel = "linear", cost = 1)
  expect_equal(
    as.character(predict(mach, array.test)@pred),
    array.test$defineCase
  )

  set.seed(12345)
  mach <- buildANN(array.train, probes = 2, size = 3, decay = 1)
  expect_equal(
    as.character(predict(mach, array.test)@pred),
    array.test$defineCase
  )
})

array.test@annot$defineCase[1] <- "Case"
array.test@annot$defineCase[10] <- "Control"

test_that("built models can detect wrong classes", {

  set.seed(12345)
  mach <- buildNB(array.train, probes = 2)
  expect_equal(
    calcStats(predict(mach, array.test), aucSkip = TRUE)$acc,
    .9
  )

  set.seed(12345)
  mach <- buildLDA(array.train, probes = 2)
  expect_equal(
    calcStats(predict(mach, array.test), aucSkip = TRUE)$acc,
    .9
  )

  set.seed(12345)
  mach <- buildSVM(array.train, probes = 2, kernel = "linear", cost = 1)
  expect_equal(
    calcStats(predict(mach, array.test), aucSkip = TRUE)$acc,
    .9
  )

  set.seed(12345)
  mach <- buildANN(array.train, probes = 2, size = 3, decay = 1)
  expect_equal(
    calcStats(predict(mach, array.test), aucSkip = TRUE)$acc,
    .9
  )
})

###########################################################
### Build and predict ExprsModule objects

set.seed(12345)

arraysMulti <- splitStratify(arrayMulti, percent.include = 80, colBy = NULL)
arrayMulti.train <- fsStats(arraysMulti[[1]], probes = 0)
arrayMulti.test <- arraysMulti[[2]]

set.seed(12345)

mach <- buildSVM(arrayMulti.train, probes = 0, kernel = "linear", cost = 1)
pred.train <- predict(mach, arrayMulti.train)
pred.test <- predict(mach, arrayMulti.test)

test_that("ExprsMulti models predict only on ExprsMulti datasets", {

  expect_error(
    predict(mach, array.train)
  )
})

mach.multi <- doMulti(arrayMulti.train, probes = 0, what = "buildSVM")

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
