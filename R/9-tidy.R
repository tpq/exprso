#' Extract Training Set
#'
#' This function extracts the training set from the result of a
#'  \code{split} method call such as \code{splitSample} or \code{splitStratify}.
#'
#' @param splitSets A two-item list. The result of a \code{split} method call.
#' @return An \code{ExprsArray} object.
#' @export
trainingSet <- function(splitSets){

  if(class(splitSets) == "list" & length(splitSets) == 2){
    return(splitSets[["array.train"]])
  }else{
    stop("Uh oh! Cannot extract a training set from this object.")
  }
}

#' Extract Validation Set
#'
#' This function extracts the validation set from the result of a
#'  \code{split} method call such as \code{splitSample} or \code{splitStratify}.
#'
#' @inheritParams trainingSet
#' @return An \code{ExprsArray} object.
#' @export
validationSet <- function(splitSets){

  if(class(splitSets) == "list" & length(splitSets) == 2){
    return(splitSets[["array.valid"]])
  }else{
    stop("Uh oh! Cannot extract a test set from this object.")
  }
}

#' @describeIn validationSet A variant of \code{validationSet}.
#' @export
testSet <- function(splitSets){

  validationSet(splitSets)
}

#' Tidy Subset Wrapper
#'
#' \code{modSubset} function provides a tidy wrapper for the \code{ExprsArray}
#'  \code{subset} method. \code{pipeSubset} provides a tidy wrapper for the
#'  \code{ExprsPipeline} \code{subset} method.
#'
#' @inheritParams arrayExprs
#' @param object An \code{ExprsArray} or \code{ExprsPipeline} object to subset.
#' @param include A character vector. Specifies which annotations in \code{colBy}
#'  to include in the subset.
#' @return An \code{ExprsArray} or \code{ExprsPipeline} object.
#' @export
modSubset <- function(object, colBy, include){

  if(inherits(object, "ExprsArray") | class(object) == "ExprsPipeline"){
    subset(object, subset = object[, colBy] %in% include)
  }else{
    stop("Uh oh! You can only use modSubset on an ExprsArray or ExprsPipeline object!")
  }
}

#' @describeIn modSubset A variant of \code{modSubset}.
#' @export
pipeSubset <- function(object, colBy, include){

  modSubset(object, colBy, include)
}
