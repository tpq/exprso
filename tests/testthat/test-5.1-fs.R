library(exprso)

set.seed(1)
fakeiris <- iris
fakeiris[,1] <- sample(fakeiris[,1])
colnames(fakeiris) <- c('bad', 'better', 'good')

binary <- exprso(fakeiris[1:100,1:3], fakeiris[1:100,5])
multi <- exprso(fakeiris[,1:3], fakeiris[,5])
cont <- exprso(fakeiris[,1:3], fakeiris[,4])

checkFS <- function(input, method, should, ...){

  args <- as.list(substitute(list(...)))[-1]
  if(identical(should, "error")){
    print(should)
    expect_error(
      do.call(method, append(args, list("object" = input)))
    )
  }else{
    print(should)
    e <- do.call(method, append(args, list("object" = input)))
    expect_equal(
      rownames(e@exprs),
      should
    )
  }
}

test_that("@preFilter order matches the @exprs order", {

  f <- fsStats(binary)
  expect_equal(
    rownames(f@exprs),
    f@preFilter[[1]]
  )

  # but now for reduction models...
  f <- fsPrcomp(binary)
  expect_false(
    identical(
      rownames(f@exprs),
      f@preFilter[[1]]
    )
  )
})

test_that("build modules work for each data type", {

  # set.seed(1);checkFS(binary, fsSample, should = c("bad", "good", "better"))
  # set.seed(1);checkFS(multi, fsSample, should = c("bad", "good", "better"))
  # set.seed(1);checkFS(cont, fsSample, should = c("bad", "good", "better"))

  checkFS(binary, fsNULL, should = c("bad", "better", "good"))
  checkFS(multi, fsNULL, should = c("bad", "better", "good"))
  checkFS(cont, fsNULL, should = c("bad", "better", "good"))

  checkFS(binary, fsANOVA, should = c("good", "better", "bad"))
  checkFS(multi, fsANOVA, should = c("good", "better", "bad"))
  checkFS(cont, fsANOVA, should = "error")

  checkFS(binary, fsInclude, should = c("bad", "good", "better"), include = c("bad", "good"))
  checkFS(multi, fsInclude, should = c("bad", "good", "better"), include = c("bad", "good"))
  checkFS(cont, fsInclude, should = c("bad", "good", "better"), include = c("bad", "good"))

  checkFS(binary, fsStats, should = c("good", "better", "bad"))
  checkFS(binary, fsStats, should = c("good", "better", "bad"), how = "ks.test")
  checkFS(binary, fsStats, should = c("good", "better", "bad"), how = "wilcox.test")
  checkFS(binary, fsStats, should = c("good", "better", "bad"), how = "var.test")
  checkFS(binary, fsStats, should = "error", how = "fake.test")
  checkFS(multi, fsStats, should = "error")
  checkFS(cont, fsStats, should = "error")

  checkFS(binary, fsCor, should = "error")
  checkFS(multi, fsCor, should = "error")
  checkFS(cont, fsCor, should = c("good", "better", "bad"))

  checkFS(binary, fsPrcomp, should = c("PC1", "PC2", "PC3"))
  checkFS(multi, fsPrcomp, should = c("PC1", "PC2", "PC3"))
  checkFS(cont, fsPrcomp, should = c("PC1", "PC2", "PC3"))

  if(requireNamespace("vegan", quietly = TRUE)){
    checkFS(binary, fsRDA, should = c("PC1", "PC2", "PC3"))
    checkFS(multi, fsRDA, should = c("PC1", "PC2", "PC3"))
    checkFS(cont, fsRDA, should = c("PC1", "PC2", "PC3"))
  }

  if(requireNamespace("limma", quietly = TRUE)){
    checkFS(binary, fsEbayes, should = c("good", "better", "bad"))
    checkFS(multi, fsEbayes, should = c("good", "better", "bad"))
    checkFS(cont, fsEbayes, should = "error")
  }

  if(requireNamespace("edgeR", quietly = TRUE)){
    checkFS(binary, fsEdger, should = c("good", "better", "bad"))
    checkFS(multi, fsEdger, should = "error")
    checkFS(cont, fsEbayes, should = "error")
  }

  if(requireNamespace("mRMRe", quietly = TRUE)){
    checkFS(binary, fsMrmre, should = c("good", "better", "bad"))
    checkFS(multi, fsMrmre, should = "error")
    checkFS(cont, fsMrmre, should = "error")
  }

  # if(requireNamespace("RankProd", quietly = TRUE)){
  #   checkFS(binary, fsRankProd, should = c("good", "better", "bad"))
  #   checkFS(multi, fsRankProd, should = "error")
  #   checkFS(cont, fsRankProd, should = "error")
  # }

  if(requireNamespace("balance", quietly = TRUE)){
    checkFS(binary, fsBalance, should = c("z1", "z2"))
    checkFS(multi, fsBalance, should = c("z1", "z2"))
    checkFS(cont, fsBalance, should = c("z1", "z2"))
  }

  checkFS(binary, fsAnnot, should = c("y", "bad", "better", "good"), colBy = "y")
  checkFS(multi, fsAnnot, should = c("y", "bad", "better", "good"), colBy = "y")
  checkFS(cont, fsAnnot, should = c("y", "bad", "better", "good"), colBy = "y")
})
