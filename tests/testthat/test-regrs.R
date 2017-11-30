library(exprso)
o <- exprso(iris[,1:3], iris[,4])
set.seed(1)
arrays <- splitSample(o)
a <- buildANN(arrays[[1]])
b <- buildRF(arrays[[1]])
c <- buildSVM(arrays[[1]])

test_that("Continuous outcome models work", {

  expect_error(
    buildLDA(o)
  )

  expect_error(
    buildNB(o)
  )

  expect_equal(
    predict(a, arrays[[2]])@actual,
    predict(b, arrays[[2]])@actual
  )

  expect_equal(
    predict(b, arrays[[2]])@actual,
    predict(c, arrays[[2]])@actual
  )
})

# ss <- ctrlSplitSet(func = "splitSample", percent.include = 67)
# fs <- ctrlFeatureSelect(func = "fsNULL", top = 0)
# gs <- ctrlGridSearch(func = "plGrid", how = "buildSVM", top = c(2, 3, 0),
#                      kernel = c("linear", "radial"))
#
# test_that("Continuous outcome pl modules work", {
#
#   set.seed(1)
#   expect_equal(
#     round(plCV(o, top = 0, how = "buildRF", fold = 2), 4),
#     .9489
#   )
#
#   set.seed(1)
#   expect_equal(
#     round(plGrid(arrays[[1]], arrays[[2]], top = c(0, 2, 3), how = "buildRF",
#                  fold = 2)@summary$valid.acc, 4),
#     c(.9645, .9244, .9643)
#   )
#
#   set.seed(1)
#   boot <- plMonteCarlo(o, B = 3, ctrlSS = ss, ctrlFS = fs, ctrlGS = gs)
#   expect_equal(
#     round(calcMonteCarlo(boot), 4),
#     .9619
#   )
#
#   set.seed(1)
#   boot <- plNested(o, fold = 2, ctrlFS = fs, ctrlGS = gs, save = FALSE)
#   expect_equal(
#     round(calcNested(boot), 4),
#     .9555
#   )
#
#   ens <- buildEnsemble(a, b, c)
#   pred <- predict(ens, o)
#   expect_equal(
#     round(calcStats(pred)$acc, 4),
#     .9286
#   )
#
#   ens <- buildEnsemble(boot)
#   pred <- predict(ens, o)
#   expect_equal(
#     round(calcStats(pred)$acc, 4),
#     .9536
#   )
# })
