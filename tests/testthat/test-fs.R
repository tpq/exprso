library(exprso)
context("fs")

load(file.path("data.RData"))

test_that("top argument to fs ExprsBinary method works", {

  expect_equal(
    rownames(fsStats(array, top = 0)@exprs),
    rownames(fsStats(array, top = c("feat4", "feat3", "feat2", "feat1"))@exprs)
  )

  expect_equal(
    rownames(fsStats(array, top = 0)@exprs),
    c("feat4", "feat1", "feat2", "feat3")
  )

  expect_equal(
    rownames(fsInclude(fsStats(array), include = "feat2")@exprs),
    c("feat2", "feat4", "feat1", "feat3")
  )

  expect_equal(
    rownames(fsStats(array, top = 0, how = "ks.test")@exprs),
    c("feat2", "feat4", "feat1", "feat3")
  )
})

test_that("@preFilter sorts features in correct order", {

  data(iris)
  y <- iris[1:100,5]
  e <- exprso(iris[1:100, 1:4], y)

  ttest <- names(
    sort(
      apply(iris[1:100,1:4], 2, function(x)
        t.test(x[y == "setosa"], x[y == "versicolor"])$p.value)
    )
  )

  expect_equal(
    fsStats(e, how = "t.test")@preFilter[[1]],
    ttest
  )

  wilcoxtest <- names(
    sort(
      apply(iris[1:100,1:4], 2, function(x)
        wilcox.test(x[y == "setosa"], x[y == "versicolor"])$p.value)
    )
  )

  expect_equal(
    fsStats(e, how = "wilcox.test")@preFilter[[1]],
    wilcoxtest
  )
})

test_that("fsANOVA matches fsStats for ExprsBinary objects", {

  expect_equal(
    fsStats(array, var.equal = TRUE)@preFilter,
    fsANOVA(array)@preFilter
  )

  expect_failure(
    expect_equal(
      fsStats(array)@preFilter,
      fsANOVA(array)@preFilter
    )
  )
})

test_that("fsANOVA works for ExprsMulti objects", {

  expect_equal(
    fsANOVA(arrayMulti)@preFilter[[1]][1],
    "feat4"
  )

  temp <- arrayMulti
  temp@exprs["feat3", temp$defineCase == 3] <- 1000
  expect_equal(
    fsANOVA(temp)@preFilter[[1]][1],
    "feat3"
  )
})
