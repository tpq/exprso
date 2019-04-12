library(exprso)
data(iris)
x <- iris[1:100,1:4]
y <- iris[1:100,5]
e <- exprso(x, y)

r1 <- fsRDA(e)
plot(r1)
m1 <- buildLR(r1, top = 1)
acc1 <- calcStats(predict(m1, e))$acc

r2 <- fsRDA(e, colBy = "defineCase")
plot(r2)
m2 <- buildLR(r2, top = 1)
acc2 <- calcStats(predict(m2, e))$acc

test_that("fsRDA will partial out colBy correctly", {

  expect_equal(
    acc1 > acc2,
    TRUE
  )
})
