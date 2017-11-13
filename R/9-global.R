#' @importFrom methods new as show
#' @importFrom plyr rbind.fill
#' @importMethodsFrom kernlab predict
#' @importMethodsFrom ROCR plot
#' @importFrom graphics plot
NULL

#' Package Check
#'
#' Checks whether the user has the required package installed.
#'  For back-end use only.
#'
#' @param package A character string. An R package.
packageCheck <- function(package){

  if(!requireNamespace(package, quietly = TRUE)){
    stop("Uh oh! This propr method depends on ", package, "!")
  }
}

#' Class Check
#'
#' Checks whether an object belongs to a specified class.
#'  For back-end use only.
#'
#' @param x An object.
#' @param what A character vector. The classes any of which \code{x} should have.
#' @param msg A string. An error message if \code{x} is not \code{what}.
classCheck <- function(x, what, msg){

  if(class(x) %in% what){ return(TRUE)
  }else if(any(sapply(what, function(cl) inherits(x, cl)))){ return(TRUE)
  }else{ stop(msg) }
}

#' Sample ExprsBinary Data
#' @usage data(array)
"array"

#' Sample ExprsMulti Data
#' @usage data(arrayMulti)
"arrayMulti"
