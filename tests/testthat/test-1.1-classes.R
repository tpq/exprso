library(exprso)

binary <- exprso(iris[1:100,1:4], iris[1:100,5])
multi <- exprso(iris[,1:4], iris[,5])
cont <- exprso(iris[,1:3], iris[,4])

pureclass <- function(data){
  expect_true(attr(class(data), "package") == "exprso")
  class(data)[1]
}

test_that("exprso objects work", {

  expect_equal(
    pureclass(binary),
    "ExprsBinary"
  )

  expect_equal(
    pureclass(multi),
    "ExprsMulti"
  )

  expect_equal(
    pureclass(cont),
    "RegrsArray"
  )

  expect_true(
    inherits(binary, "ExprsArray")
  )

  expect_true(
    inherits(multi, "ExprsArray")
  )

  expect_true(
    inherits(cont, "ExprsArray")
  )

  expect_equal(
    pureclass(buildRF(binary)),
    "ExprsMachine"
  )

  expect_equal(
    pureclass(buildRF(multi)),
    "ExprsModule"
  )

  expect_equal(
    pureclass(buildRF(cont)),
    "RegrsModel"
  )

  expect_true(
    inherits(buildRF(binary), "ExprsModel")
  )

  expect_true(
    inherits(buildRF(multi), "ExprsModel")
  )

  expect_true(
    inherits(buildRF(cont), "ExprsModel")
  )

  expect_equal(
    pureclass(predict(buildRF(binary), binary)),
    "ExprsPredict"
  )

  expect_equal(
    pureclass(predict(buildRF(multi), multi)),
    "MultiPredict"
  )

  expect_equal(
    pureclass(predict(buildRF(cont), cont)),
    "RegrsPredict"
  )
})
