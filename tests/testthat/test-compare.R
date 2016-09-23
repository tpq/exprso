library(exprso)
context("compare")

###########################################################
### Compare ExprsBinary and ExprsMulti objects

load(file.path("data.RData"))

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
