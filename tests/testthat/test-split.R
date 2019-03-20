library(exprso)
suppressWarnings(RNGversion("3.5.0"))

###########################################################
### Test ExprsBinary and ExprsMulti imports

set.seed(1235)

df.a <- data.frame(
  "id" = 1:10,
  "class" = rep("a", 10),
  "sex" = c(rep("M", 5), rep("F", 5)),
  "feat1" = rnorm(10, mean = 10, sd = 1),
  "feat2" = rnorm(10, mean = 20, sd = 5),
  "feat3" = rnorm(10, mean = 5, sd = 1)
)

df.b <- data.frame(
  "id" = 11:30,
  "class" = rep("b", 20),
  "sex" = c(rep("M", 10), rep("F", 10)),
  "feat1" = rnorm(20, mean = 20, sd = 5),
  "feat2" = rnorm(20, mean = 10, sd = 1),
  "feat3" = rnorm(20, mean = 5, sd = 1)
)

df.c <- data.frame(
  "id" = 31:40,
  "class" = rep("c", 10),
  "sex" = c(rep("M", 3), rep("F", 7)),
  "feat1" = rnorm(10, mean = 15, sd = 3),
  "feat2" = rnorm(10, mean = 15, sd = 3),
  "feat3" = rnorm(10, mean = 5, sd = 1)
)

df <- do.call(rbind, list(df.a, df.b, df.c))

tempFile <- tempfile()
write.table(df, file = tempFile, sep = "\t")

array <-
  arrayExprs(tempFile, begin = 4, colID = "id", colBy = "class",
             include = list("a", "b"))

arrayMulti <-
  arrayExprs(tempFile, begin = 4, colID = "id", colBy = "class",
             include = list("a", "b", "c"))

test_that("ExprsBinary imports correctly", {

  expect_equal(
    as.character(class(array)),
    "ExprsBinary"
  )

  expect_equal(
    unique(array[array$defineCase == "Control", "class"]),
    "a"
  )

  expect_equal(
    unique(array[array$defineCase == "Case", "class"]),
    "b"
  )
})

test_that("ExprsMulti imports correctly", {

  expect_equal(
    as.character(class(arrayMulti)),
    "ExprsMulti"
  )

  expect_equal(
    unique(arrayMulti[arrayMulti$defineCase == 1, "class"]),
    "a"
  )

  expect_equal(
    unique(arrayMulti[arrayMulti$defineCase == 2, "class"]),
    "b"
  )

  expect_equal(
    unique(arrayMulti[arrayMulti$defineCase == 3, "class"]),
    "c"
  )
})

###########################################################
### Test splitStratify without colBy

arrays <-
  splitStratify(array, percent.include = 50, colBy = NULL)

arraysMulti <-
  splitStratify(arrayMulti, percent.include = 50, colBy = NULL)

test_that("splitStratify correctly splits ExprsBinary objects", {

  expect_equal(
    nrow(arrays[[1]]@annot),
    10
  )

  expect_equal(
    sum(arrays[[2]]$defineCase == "Control"),
    5
  )

  expect_equal(
    sum(arrays[[2]]$defineCase == "Case"),
    15
  )
})

test_that("splitStratify correctly splits ExprsMulti objects", {

  expect_equal(
    sum(arraysMulti[[1]]$defineCase == 1),
    5
  )

  expect_equal(
    sum(arraysMulti[[2]]$defineCase == 2),
    15
  )

  expect_equal(
    sum(arraysMulti[[2]]$defineCase == 3),
    5
  )
})

###########################################################
### Test splitStratify with colBy

arrays <-
  splitStratify(array, percent.include = 100, colBy = c("sex"), bin = c(FALSE))

arraysMulti <-
  splitStratify(arrayMulti, percent.include = 100, colBy = c("sex"), bin = c(FALSE))

test_that("splitStratify correctly splits ExprsBinary objects with colBy", {

  expect_equal(
    nrow(arrays[[1]]@annot),
    20
  )

  expect_equal(
    nrow(arrays[[2]]@annot),
    10
  )

  expect_equal(
    sum(arrays[[1]]$defineCase == "Control" & arrays[[1]]$sex == "M"),
    sum(arrays[[1]]$defineCase == "Case" & arrays[[1]]$sex == "M")
  )

  expect_equal(
    sum(arrays[[1]]$defineCase == "Control" & arrays[[1]]$sex == "F"),
    sum(arrays[[1]]$defineCase == "Case" & arrays[[1]]$sex == "F")
  )
})

test_that("splitStratify correctly splits ExprsMulti objects with colBy", {

  expect_equal(
    nrow(arraysMulti[[1]]@annot),
    24
  )

  expect_equal(
    nrow(arraysMulti[[2]]@annot),
    16
  )

  expect_equal(
    sum(arraysMulti[[1]]$defineCase == 1 & arraysMulti[[1]]$sex == "M"),
    3
  )

  expect_equal(
    sum(arraysMulti[[1]]$defineCase == 2 & arraysMulti[[1]]$sex == "M"),
    3
  )

  expect_equal(
    sum(arraysMulti[[1]]$defineCase == 3 & arraysMulti[[1]]$sex == "M"),
    3
  )

  expect_equal(
    sum(arraysMulti[[1]]$defineCase == 1 & arraysMulti[[1]]$sex == "F"),
    5
  )

  expect_equal(
    sum(arraysMulti[[1]]$defineCase == 2 & arraysMulti[[1]]$sex == "F"),
    5
  )

  expect_equal(
    sum(arraysMulti[[1]]$defineCase == 3 & arraysMulti[[1]]$sex == "F"),
    5
  )
})
