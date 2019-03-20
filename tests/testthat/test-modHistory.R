library(exprso)
data(iris)
e <- exprso(iris[1:100,1:4], iris[1:100,5])

A <- fsPCA(e)
B <- fsRDA(e)
C <- fsBalance(e)

test_that("modHistory reconstructs reduced dimensions", {

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
})
