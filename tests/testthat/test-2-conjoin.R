library(exprso)

binary <- exprso(iris[1:100,1:4], iris[1:100,5])
multi <- exprso(iris[,1:4], iris[,5])
cont <- exprso(iris[,1:3], iris[,4])

build.binary <- plGrid(binary, how = "buildSVM", cost = 1:7, top = 0)
build.multi <- plGrid(binary, how = "buildSVM", cost = 1:7, top = 0)
build.cont <- plGrid(binary, how = "buildSVM", cost = 1:7, top = 0)

checkConjoin <- function(object){

  a <- object[1:5,]
  b <- object[6:7,]

  ab <- conjoin(a, b)
  ba <- conjoin(b, a)
  if(class(object) == "ExprsPipeline"){
    ab@summary <- ab@summary[,-1] # get rid of boot column
    ba@summary <- ba@summary[,-1]
  }

  expect_equal(
    ab,
    object[1:7,]
  )

  expect_equal(
    ba,
    object[c(6:7, 1:5),]
  )
}

test_that("conjoin method works correctly", {

  checkConjoin(binary)
  checkConjoin(multi)
  checkConjoin(cont)

  checkConjoin(build.binary)
  checkConjoin(build.multi)
  checkConjoin(build.cont)

  set.seed(1)
  c.mac <- conjoin(buildSVM(binary), buildLASSO(binary),
                   buildRF(binary), buildDT(binary))
  set.seed(1)
  c.ens <- conjoin(conjoin(buildSVM(binary), buildLASSO(binary)),
                   conjoin(buildRF(binary), buildDT(binary)))
  expect_equal(
    c.mac,
    c.ens
  )
})
