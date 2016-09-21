library(exprso)
library(magrittr)
context("mod")

###########################################################
### Test modSwap, modCluster, and modSubset

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
  "sex" = c(rep("M", 4), rep("F", 6)),
  "feat1" = rnorm(10, mean = 15, sd = 3),
  "feat2" = rnorm(10, mean = 15, sd = 3),
  "feat3" = rnorm(10, mean = 5, sd = 1),
  "feat4" = rnorm(10, mean = 30, sd = 1)
)

df <- do.call(rbind, list(df.a, df.b, df.c))

tempFile <- tempfile()
write.table(df, file = tempFile, sep = "\t")

array <-
  arrayExprs(tempFile, probes.begin = 4, colID = "id", colBy = "class",
             include = list("a", "b"))

arrays <- splitStratify(array, percent.include = 50, colBy = "sex")
array.train <- arrays[[1]]
spill <- splitStratify(arrays[[2]], percent.include = 67, colBy = "sex")
array.test <- spill[[1]]
spillover <- spill[[2]]

test_that("modSwap and subset work without any error", {

  set.seed(1)
  expect_equal(
    array.train %>% modSwap(percent = 50) %>% subset(select = "mutated") %>% '$'("mutated"),
    c(0, 1, 0, 1, 0, 0, 0, 0)
  )

  set.seed(1)
  expect_equal(

    array.train %>% modSwap(how = spillover, percent = 50) %>% '$'("mutated"),
    c(0, 1, 0, 1, 0, 0, 0, 0)
  )
})

test_that("modCluster and modSubset work without any error", {

  expect_equal(

    array.train %>% modCluster(how = "hclust") %>% modSubset(colBy = "cluster", include = 1),
    array.train %>% modCluster(how = "kmeans") %>% modSubset(colBy = "cluster", include = 1)
  )

  expect_equal(

    array.train %>% modCluster(how = "agnes", k = 3) %>% modSubset(colBy = "cluster", include = 3),
    array.train %>% modCluster(how = "clara", k = 3) %>% modSubset(colBy = "cluster", include = 3)
  )

  expect_equal(

    array.train %>% modCluster(how = "diana"),
    array.train %>% modCluster(how = "fanny")
  )

  expect_equal(

    array.train %>% modCluster(how = "hclust"),
    array.train %>% modCluster(how = "pam")
  )
})
