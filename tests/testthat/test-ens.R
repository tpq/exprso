library(exprso)
context("ens")

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
  arrayExprs(tempFile, begin = 4, colID = "id", colBy = "class",
             include = list("a", "b"))

###########################################################
### Check plMonteCarlo

array@annot$defineCase[1:2] <- "Case"
array@annot$defineCase[29:30] <- "Control"

set.seed(12345)
arrays <- splitStratify(array, percent.include = 50, colBy = NULL)

# Build ensemble with plGrid
set.seed(12345)
pl <- plGrid(arrays[[1]], top = 0, how = "buildSVM", fold = NULL,
             kernel = "linear", cost = 10^(c(-3, 1, 3)))
ens1 <- buildEnsemble(pl)

# Build ensemble manually
set.seed(12345)
mach1 <- buildSVM(arrays[[1]], top = 0, kernel = "linear", cost = 10^-3)
mach2 <- buildSVM(arrays[[1]], top = 0, kernel = "linear", cost = 10^1)
mach3 <- buildSVM(arrays[[1]], top = 0, kernel = "linear", cost = 10^3)
ens2 <- buildEnsemble(mach1, mach2, mach3, "notMachine")

test_that("buildEnsemble and conjoin work the same", {

  expect_equal(
    conjoin(mach1, mach2, mach3, "NOTmachine"),
    ens2
  )
})

test_that("buildEnsemble methods for pl and model match", {

  expect_equal(
    predict(ens1, arrays[[2]]),
    predict(ens2, arrays[[2]])
  )
})

ens3 <- buildEnsemble(pl, colBy = "train.acc", how = .75)
ens4 <- buildEnsemble(mach3, mach2)

test_that("buildEnsemble calls pipeFilter correctly", {

  expect_equal(
    predict(ens3, arrays[[2]]),
    predict(ens4, arrays[[2]])
  )
})
