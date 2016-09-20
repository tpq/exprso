library(exprso)
context("plGrid")
context("plCV")

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

###########################################################
### Check plGrid

arrays <- splitStratify(array, percent.include = 50, colBy = NULL)
array.train <- arrays[[1]]
array.test <- arrays[[2]]

pl <- plGrid(array.train,
             array.test,
             probes = 0,
             how = "buildLDA",
             fold = NULL,
             method = c("mle",
                        "mve")
)

test_that("individual @machs match expectations", {

  mach <- buildLDA(array.train, probes = 0, method = "mle")
  expect_equal(
    predict(pl@machs[[1]], array.test),
    predict(mach, array.test)
  )

  mach <- buildLDA(array.train, probes = 0, method = "mve")
  expect_equal(
    predict(pl@machs[[2]], array.test),
    predict(mach, array.test)
  )
})

test_that("@summary match expectations", {

  mach <- buildLDA(array.train, probes = 0, method = "mle")
  expect_equal(
    matrix(pl@summary[1, c("valid.acc", "valid.sens", "valid.spec", "valid.auc")]),
    matrix(calcStats(predict(mach, array.test)))
  )

  mach <- buildLDA(array.train, probes = 0, method = "mve")
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

acc <- plCV(array.train, probes = 0, how = "buildLDA", fold = 10, method = "mle")

array.train@annot$defineCase[c(1:2)] <- "Case"
array.train@annot$defineCase[c(19:20)] <- "Control"

acc.off <- plCV(array.train, probes = 0, how = "buildLDA", fold = 0, method = "mle")
pl <- plGrid(array.train,
             array.test,
             probes = 0,
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
