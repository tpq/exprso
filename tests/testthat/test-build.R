library(exprso)

binary <- exprso(iris[1:100,1:4], iris[1:100,5])
multi <- exprso(iris[,1:4], iris[,5])
cont <- exprso(iris[,1:3], iris[,4])

checkBuild <- function(input, classifier, should){

  if(should == "error"){
    print(should)
    expect_error(
      do.call(classifier, list("object" = input))
    )
  }else{
    print(should)
    m <- do.call(classifier, list("object" = input))
    expect_equal(
      round(calcStats(predict(m, input))$acc, 3),
      round(should, 3)
    )
  }
}

test_that("build modules work for each data type", {

  # NAIVE BAYES
  set.seed(1)
  checkBuild(binary, buildNB, should = 1)
  checkBuild(multi, buildNB, should = .96)
  checkBuild(cont, buildNB, should = "error")

  # LINEAR DISCRIMINANT ANALYSIS
  set.seed(1)
  checkBuild(binary, buildLDA, should = 1)
  checkBuild(multi, buildLDA, should = .98)
  checkBuild(cont, buildLDA, should = "error")

  # SUPPORT VECTOR MACHINE
  set.seed(1)
  checkBuild(binary, buildSVM, should = 1)
  checkBuild(multi, buildSVM, should = .9667)
  checkBuild(cont, buildSVM, should = .9378)

  # LM / GLM / LR
  set.seed(1)
  checkBuild(binary, buildLM, should = "error")
  checkBuild(multi, buildLM, should = "error")
  checkBuild(cont, buildLM, should = .9379)
  set.seed(1)
  checkBuild(binary, buildGLM, should = 1)
  checkBuild(multi, buildGLM, should = "error")
  checkBuild(cont, buildGLM, should = .9379)

  # LASSO
  set.seed(1)
  checkBuild(binary, buildLASSO, should = 1)
  checkBuild(multi, buildLASSO, should = .953)
  checkBuild(cont, buildLASSO, should = .6690)

  # NEURAL NETS
  set.seed(1)
  checkBuild(binary, buildANN, should = 1)
  checkBuild(multi, buildANN, should = .6667)
  #checkBuild(cont, buildANN, should = NA)

  # DECISION TREES
  set.seed(1)
  checkBuild(binary, buildDT, should = 1)
  checkBuild(multi, buildDT, should = .96)
  checkBuild(cont, buildDT, should = .9336)

  # RANDOM FORESTS
  set.seed(1)
  checkBuild(binary, buildRF, should = 1)
  checkBuild(multi, buildRF, should = 1)
  checkBuild(cont, buildRF, should = .9778)

  # FRB
  set.seed(1)
  checkBuild(binary, buildFRB, should = 1)
  checkBuild(multi, buildFRB, should = .953)
  checkBuild(cont, buildFRB, should = .932)
})
