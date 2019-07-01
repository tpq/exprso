library(exprso)
data(iris)
e <- exprso(iris[1:100,1:4], iris[1:100,5])

A <- fsPCA(e)
B <- fsRDA(e)
C <- fsBalance(e)
D <- fsAnnot(e, colBy = "y")

test_that("modHistory reconstructs reduced dimensions", {

  expect_equal(
    fsPrcomp(e),
    fsPCA(e)
  )

  expect_equal(
    A@exprs,
    modHistory(e, A)@exprs
  )

  expect_equal(
    B@exprs,
    modHistory(e, B)@exprs
  )

  expect_equal(
    C@exprs,
    modHistory(e, C)@exprs
  )

  expect_equal(
    D@exprs,
    modHistory(e, D)@exprs
  )
})
