library(exprso)
context("fs")

###########################################################
### Prepare ExprsBinary and ExprsMulti objects

set.seed(1235)

df.a <- data.frame(
  "id" = 1:10,
  "class" = rep("a", 10),
  "sex" = c(rep("M", 5), rep("F", 5)),
  "feat1" = rnorm(10, mean = 10, sd = 1),
  "feat2" = rnorm(10, mean = 20, sd = 5),
  "feat3" = rnorm(10, mean = 5, sd = 1),
  "feat4" = rnorm(10, mean = 40, sd = 1)
)

df.b <- data.frame(
  "id" = 11:30,
  "class" = rep("b", 20),
  "sex" = c(rep("M", 10), rep("F", 10)),
  "feat1" = rnorm(20, mean = 20, sd = 5),
  "feat2" = rnorm(20, mean = 10, sd = 1),
  "feat3" = rnorm(20, mean = 5, sd = 1),
  "feat4" = rnorm(20, mean = 20, sd = 1)
)

df.c <- data.frame(
  "id" = 31:40,
  "class" = rep("c", 10),
  "sex" = c(rep("M", 3), rep("F", 7)),
  "feat1" = rnorm(10, mean = 15, sd = 3),
  "feat2" = rnorm(10, mean = 15, sd = 3),
  "feat3" = rnorm(10, mean = 5, sd = 1),
  "feat4" = rnorm(10, mean = 30, sd = 1)
)

df <- do.call(rbind, list(df.a, df.b, df.c))

tempFile <- tempfile()
write.table(df, file = tempFile, sep = "\t")

array <-
  arrayRead(tempFile, probes.begin = 4, colID = "id", colBy = "class",
            include = list("a", "b"))

arrayMulti <-
  arrayRead(tempFile, probes.begin = 4, colID = "id", colBy = "class",
            include = list("a", "b", "c"))

###########################################################
### Test ExprsBinary feature selection

test_that("probes argument to fs ExprsBinary method works", {

  expect_equal(
    rownames(fsStats(array, probes = 0)@exprs),
    rownames(fsStats(array, probes = c("feat4", "feat3", "feat2", "feat1"))@exprs)
  )

  expect_equal(
    rownames(fsStats(array, probes = 0)@exprs),
    c("feat4", "feat1", "feat2", "feat3")
  )

  expect_equal(
    rownames(fsStats(array, probes = 0, how = "ks.test")@exprs),
    rownames(fsStats(array, probes = c("feat4", "feat3", "feat2", "feat1"), how = "ks.test")@exprs)
  )

  expect_equal(
    rownames(fsStats(array, probes = 0, how = "ks.test")@exprs),
    c("feat2", "feat4", "feat1", "feat3")
  )

  expect_equal(
    rownames(fsEbayes(array, probes = 0)@exprs),
    rownames(fsEbayes(array, probes = c("feat4", "feat3", "feat2", "feat1"))@exprs)
  )

  expect_equal(
    rownames(fsEbayes(array, probes = 0)@exprs),
    c("feat1", "feat3", "feat2", "feat4")
  )

  expect_equal(
    rownames(fsMrmre(array, probes = 0)@exprs),
    rownames(fsMrmre(array, probes = c("feat4", "feat3", "feat2", "feat1"))@exprs)
  )

  expect_equal(
    rownames(fsMrmre(array, probes = 0)@exprs),
    c("feat4", "feat2", "feat1", "feat3")
  )

  expect_equal(
    rownames(fsPathClassRFE(array, probes = 0)@exprs),
    c("feat2", "feat4")
  )
})

###########################################################
### Test ExprsMulti feature selection

fsStats.multi <- doMulti(arrayMulti, probes = 0, what = "fsStats")

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

test_that("probes argument to fs ExprsBinary method works", {

  expect_equal(
    rownames(fsStats(arrayMulti, probes = 0)@exprs),
    c("feat4", "feat1", "feat2", "feat3")
  )

  expect_equal(

    rownames(fsStats(arrayMulti, probes = 0)@exprs),
    rownames(fsStats(arrayMulti, probes = c("feat4", "feat3", "feat2", "feat1"))@exprs)
  )
})
