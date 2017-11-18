#' @importFrom methods new as show
#' @importFrom stats aov t.test ks.test prcomp
#' @importFrom utils write.csv
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
    stop("Uh oh! This propr method depends on ", package, ".")
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

#' Build an args List
#' @param ... Arguments passed down from a calling function.
getArgs <- function(...){

  args <- as.list(substitute(list(...)))[-1]
  return(args)
}

#' Set an args List Element to Default Value
#' @param what The name of the argument.
#' @param as The value to set it as.
#' @param args An args list. The result of \code{\link{getArgs}}.
defaultArg <- function(what, as, args){

  if(!what %in% names(args)){
    cat("Setting", what, "to", as.character(as), "(default behavior, override explicitly)...\n")
    as <- list(as); names(as) <- what
    args <- append(args, as)
  }

  return(args)
}

#' Force an args List Element to Value
#' @inheritParams defaultArg
forceArg <- function(what, as, args){

  if(!what %in% names(args)){
    cat("Setting", what, "to", as.character(as), "(forced behavior, cannot override)...\n")
    as <- list(as); names(as) <- what
    args <- append(args, as)
  }else{
    if(args[[what]] == as){
      cat(paste0("Uh oh! This function requires ", what, " = ", as,
                 ". Setting ", what, " to ", as, "...\n"))
      args[[what]] <- as
    }
  }

  return(args)
}

#' Sample ExprsBinary Data
#' @usage data(array)
"array"

#' Sample ExprsMulti Data
#' @usage data(arrayMulti)
"arrayMulti"
