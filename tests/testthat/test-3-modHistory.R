library(exprso)
data(iris)
e <- exprso(iris[1:100,1:4], iris[1:100,5])

test_that("modHistory reconstructs reduced dimensions", {

  expect_equal(
    fsPrcomp(e),
    fsPCA(e)
  )

  A <- fsPCA(e)
  expect_equal(
    A@exprs,
    modHistory(e, A)@exprs
  )

  B <- fsRDA(e)
  expect_equal(
    B@exprs,
    modHistory(e, B)@exprs
  )

  C <- fsBalance(e)
  expect_equal(
    C@exprs,
    modHistory(e, C)@exprs
  )

  D <- fsAnnot(e, colBy = "y")
  expect_equal(
    D@exprs,
    modHistory(e, D)@exprs
  )

  if(requireNamespace("amalgam", quietly = TRUE)){

    E <- fsAmalgam(e)
    expect_equal(
      E@exprs,
      modHistory(e, E)@exprs
    )

    F <- fsAmalgam(e, n.amalgams = 4, asSLR = TRUE)
    expect_equal(
      F@exprs,
      modHistory(e, F)@exprs
    )
  }

  G <- fsPRA(e, nclust = 3)
  expect_equal(
    G@exprs,
    modHistory(e, G)@exprs
  )
})
