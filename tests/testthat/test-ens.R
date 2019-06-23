library(exprso)

binary <- exprso(iris[1:100,1:4], iris[1:100,5])
multi <- exprso(iris[,1:4], iris[,5])
cont <- exprso(iris[,1:3], iris[,4])

checkEnsemble <- function(ens, data, should){

  print(should)
  expect_equal(
    round(calcStats(predict(ens, data))$acc, 3),
    round(should, 3)
  )
}

test_that("buildEnsemble is actually the same as conjoin", {

  set.seed(1); b <- buildEnsemble(buildSVM(binary), buildLASSO(binary), buildRF(binary))
  set.seed(1); j <- conjoin(buildSVM(binary), buildLASSO(binary), buildRF(binary))
  expect_equal(b, j)
})

test_that("buildEnsemble from argument works for each data type", {

  set.seed(1)
  ens.binary <- buildEnsemble(buildSVM(binary), buildLASSO(binary), buildRF(binary))
  checkEnsemble(ens.binary, binary, should = 1)

  set.seed(1)
  ens.multi <- buildEnsemble(buildSVM(multi), buildLASSO(multi), buildRF(multi))
  checkEnsemble(ens.multi, multi, should = .9733)

  set.seed(1)
  ens.cont <- buildEnsemble(buildSVM(cont), buildLASSO(cont), buildRF(cont))
  checkEnsemble(ens.cont, cont, should = .9649)
})

test_that("buildEnsemble from pl works for each data type", {

  set.seed(1)
  ens.binary <- buildEnsemble(plGrid(binary, how = "buildSVM", top = c(1, 2, 3)))
  checkEnsemble(ens.binary, binary, should = 1)

  set.seed(1)
  ens.multi <- buildEnsemble(plGrid(multi, how = "buildSVM", top = c(1, 2, 3)))
  checkEnsemble(ens.multi, multi, should = .8333)

  set.seed(1)
  ens.cont <- buildEnsemble(plGrid(cont, how = "buildSVM", top = c(1, 2, 3)))
  checkEnsemble(ens.cont, cont, should = .8472)
})
