library(exprso)
context("compare")

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
  "sex" = c(rep("M", 4), rep("F", 6)),
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
### Compare ExprsBinary and ExprsMulti objects

arrays <- splitStratify(array, percent.include = 50, colBy = "sex")
array.train <- arrays[[1]]
array.test <- splitStratify(arrays[[2]], percent.include = 100, colBy = "sex")[[1]]

test_that("compare works for ExprsBinary objects", {

  expect_equal(
    compare(array.train)[[1]],
    compare(array.train, array.test)[[1]]
  )

  expect_equal(
    compare(array.train)[[1]],
    c("id" = TRUE, "class" = TRUE, "sex" = FALSE)
  )

  expect_equal(
    compare(array.train, array.test)[[3]],
    c("id" = FALSE, "class" = FALSE, "sex" = FALSE, "defineCase" = FALSE)
  )
})

arrays <- splitStratify(arrayMulti, percent.include = 50, colBy = "sex")
array.train <- arrays[[1]]
array.test <- splitStratify(arrays[[2]], percent.include = 100, colBy = "sex")[[1]]

test_that("compare works for ExprsMulti objects", {

  expect_equal(
    compare(array.train)[[1]],
    compare(array.train, array.test)[[1]]
  )

  expect_equal(
    compare(array.train)[[1]],
    c("id" = TRUE, "class" = TRUE, "sex" = FALSE)
  )

  expect_equal(
    compare(array.train, array.test)[[3]],
    c("id" = FALSE, "class" = FALSE, "sex" = FALSE, "defineCase" = FALSE)
  )
})
