library(exprso)
context("fs")

if(requireNamespace("limma", quietly = TRUE) &
   requireNamespace("edgeR", quietly = TRUE) &
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

obj1 <- exprso(iris[1:100, 1:4], iris[1:100, 5])
obj2 <- exprso(iris[,1:4], iris[,5])
obj3 <- exprso(iris[,1:3], iris[,4])

test_that("fsANOVA works", {

  expect_equal(
    fsANOVA(obj1)@preFilter[[1]],
    fsStats(obj1, var.equal = TRUE)@preFilter[[1]]
  )

  expect_equal(
    fsANOVA(obj1)@preFilter[[1]],
    fsANOVA(obj2)@preFilter[[1]]
  )

  expect_error(
    fsANOVA(obj3)
  )
})

test_that("fsCor works", {

  expect_error(
    fsCor(obj1)
  )

  expect_error(
    fsCor(obj2)
  )

  expect_equal(
    fsCor(obj3)@preFilter[[1]],
    c("Petal.Length", "Sepal.Length", "Sepal.Width")
  )
})
