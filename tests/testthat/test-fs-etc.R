library(exprso)
context("fs")

###########################################################
### Test ExprsBinary feature selection

if(requireNamespace("limma", quietly = TRUE) &
   requireNamespace("mRMRe", quietly = TRUE) &
   requireNamespace("pathClass", quietly = TRUE) &
   requireNamespace("propr", quietly = TRUE)){

  load(file.path("data.RData"))

  test_that("top argument to fs ExprsBinary method works", {

    expect_equal(
      rownames(fsEbayes(array, top = 0)@exprs),
      rownames(fsEbayes(array, top = c("feat4", "feat3", "feat2", "feat1"))@exprs)
    )

    expect_equal(
      rownames(fsEbayes(array, top = 0)@exprs),
      c("feat4", "feat2", "feat3", "feat1")
    )

    expect_equal(
      rownames(fsEdger(array, top = 0)@exprs),
      c("feat1", "feat4", "feat2", "feat3")
    )

    expect_equal(
      rownames(fsMrmre(array, top = 0)@exprs),
      rownames(fsMrmre(array, top = c("feat4", "feat3", "feat2", "feat1"))@exprs)
    )

    expect_equal(
      rownames(fsMrmre(array, top = 0)@exprs),
      c("feat4", "feat2", "feat1", "feat3")
    )

    expect_equal(
      rownames(fsPathClassRFE(array, top = 0)@exprs),
      c("feat2", "feat4")
    )

    expect_equal(
      rownames(fsPropd(array, top = 0)@exprs),
      c("feat1", "feat4", "feat2", "feat3")
    )
  })
}
