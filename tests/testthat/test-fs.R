library(exprso)
context("fs")

###########################################################
### Test ExprsBinary feature selection

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
    rownames(fsStats(array, top = 0, how = "ks.test")@exprs),
    c("feat2", "feat4", "feat1", "feat3")
  )

  expect_equal(
    rownames(fsEbayes(array, top = 0)@exprs),
    rownames(fsEbayes(array, top = c("feat4", "feat3", "feat2", "feat1"))@exprs)
  )

  expect_equal(
    rownames(fsEbayes(array, top = 0)@exprs),
    c("feat1", "feat3", "feat2", "feat4")
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

###########################################################
### Test ExprsMulti feature selection

fsStats.multi <- doMulti(arrayMulti, top = 0, method = "fsStats")

test_that("doMulti performs 1 vs. all fs", {

  expect_equal(
    rownames(fsStats.multi[[1]]@exprs)[1],
    "feat4"
  )

  expect_equal(
    rownames(fsStats.multi[[1]]@exprs)[4],
    "feat3"
  )

  expect_equal(
    rownames(fsStats.multi[[2]]@exprs)[1],
    "feat4"
  )

  expect_equal(
    rownames(fsStats.multi[[2]]@exprs)[4],
    "feat3"
  )

  expect_equal(
    rownames(fsStats.multi[[3]]@exprs)[1],
    "feat4"
  )

  expect_equal(
    rownames(fsStats.multi[[3]]@exprs)[4],
    "feat3"
  )
})

test_that("reRank does not yield egregious error", {

  expect_equal(
    reRank(fsStats.multi)[1],
    "feat4"
  )

  expect_equal(
    reRank(fsStats.multi)[4],
    "feat3"
  )
})
