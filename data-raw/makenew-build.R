###
# Set up data each with binary, mulit-class, and continuous outcomes
###

library(exprso)
data(iris)
e1 <- exprso(iris[1:100, 1:4], iris[1:100, 5])
e1
e2 <- e <- exprso(iris[, 1:4], iris[, 5])
e2
e3 <- exprso(iris[, 1:3], iris[, 4])
e3

###
# Bring workhorses into active environment
###

fs. <- exprso:::fs.
build. <- exprso:::build.
getArgs <- exprso:::getArgs
defaultArg <- exprso:::defaultArg
classCheck <- exprso:::classCheck
forceArg <- exprso:::forceArg

###
# Create new module here
###

#' Build Logistic Regression Model
#'
#' \code{buildLR} builds a model using the \code{glm} function.
#'
#' @inheritParams build.
#' @return Returns an \code{ExprsModel} object.
#' @export
buildLR <- function(object, top = 0, ...){ # args to glm

  classCheck(object, c("ExprsBinary", "ExprsMulti"),
             "This build method only works for classification tasks.")

  build.(object, top,
         uniqueFx = function(data, labels, ...){

           # Perform GLM via ~ method
           args <- getArgs(...)
           args <- forceArg("family", "binomial", args)
           df <- data.frame(data, "defineCase" = as.numeric(labels) - 1)
           args <- append(list("formula" = defineCase ~ ., "data" = df), args)
           do.call(stats::glm, args)
         }, ...)
}

###
# Test new module here (pre-implementation)
###

m1 <- buildLR(e1)
m2 <- buildLR(e2)
m3 <- buildLR(e3)

predict(m1@mach, as.data.frame(t(e1@exprs)))
predict(m3@mach, as.data.frame(t(e3@exprs)))

class(m1@mach)
class(m3@mach)

###
# Test new module here (post-implementation)
###

library(exprso)
m1 <- buildLR(e1)
m2 <- buildLR(e2)
m3 <- buildLR(e3)
predict(m1, e1)
predict(m2, e2)
predict(m3, e3)
