library(exprso)
o <- exprso(iris[,1:3], iris[,4])
arrays <- splitSample(o)
a <- buildANN(arrays[[1]])
b <- buildRF(arrays[[1]])
c <- buildSVM(arrays[[1]])

test_that("Continuous outcome models work", {

  expect_error(
    buildLDA(array)
  )

  expect_error(
    buildNB(array)
  )

  expect_equal(
    predict(a, arrays[[2]])@actual,
    predict(b, arrays[[2]])@actual
  )

  expect_equal(
    predict(b, arrays[[2]])@actual,
    predict(c, arrays[[2]])@actual
  )
})
