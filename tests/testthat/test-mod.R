library(exprso)
library(magrittr)
context("mod")

###########################################################
### Check modSwap, modSubset, and modCluster

load(file.path("data.RData"))

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

    array.train %>% modCluster(how = "diana"),
    array.train %>% modCluster(how = "fanny")
  )

  expect_equal(

    array.train %>% modCluster(how = "hclust"),
    array.train %>% modCluster(how = "pam")
  )
})
