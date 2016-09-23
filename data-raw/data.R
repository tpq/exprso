set.seed(1235) # changing seed may break tests!

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

arrayMulti <-
  arrayExprs(tempFile, begin = 4, colID = "id", colBy = "class",
             include = list("a", "b", "c"))

devtools::use_data(array, arrayMulti)
