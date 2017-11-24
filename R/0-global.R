#' @importFrom methods new as show
#' @importFrom stats aov t.test ks.test prcomp cor
#' @importFrom stats median quantile
#' @importFrom utils write.csv
#' @importFrom plyr rbind.fill
#' @importMethodsFrom kernlab predict
#' @importMethodsFrom ROCR plot
#' @importFrom graphics plot
NULL

#' Sample ExprsBinary Data
#' @usage data(array)
"array"

#' Sample ExprsMulti Data
#' @usage data(arrayMulti)
"arrayMulti"

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

#' Build Argument Grid
#'
#' This function builds an argument grid from any number of arguments.
#'  Used to prepare a grid-search for the \code{plGrid} and
#'  \code{plGridMulti} functions.
#'
#' @param array.train The \code{array.train} argument as fed to \code{plGrid}.
#' @param top The \code{top} argument as fed to \code{plGrid}.
#' @param how The \code{how} argument as fed to \code{plGrid}.
#' @param ... Additional arguments as fed to \code{plGrid}.
makeGridFromArgs <- function(array.train, top, how, ...){

  if(is.numeric(top)){
    if(any(top > nrow(array.train@exprs))){
      message("At least one 'top' index is too large. Using all features instead.")
      top[top > nrow(array.train@exprs)] <- nrow(array.train@exprs)
    }
    top <- unique(top)
  }else if(is.character(top)){
    top <- as.list(top)
  } # if list, do nothing here

  args <- getArgs(...)
  args <- append(list("top" = top), lapply(args, eval))
  grid <- expand.grid(args, stringsAsFactors = FALSE)

  if(how == "buildSVM"){
    if(!"kernel" %in% names(args)) grid$kernel <- "linear"
    if(!"cost" %in% names(args)) grid$cost <- 1
    if("radial" %in% grid$kernel){
      if(!"gamma" %in% names(args)) grid$gamma <- 0.1
      grid[grid$kernel %in% "linear", "gamma"] <- NA
    }
    if("polynomial" %in% grid$kernel){
      if(!"degree" %in% names(args)) grid$degree <- 2
      if(!"coef0" %in% names(args)) grid$coef0 <- 1
      grid[grid$kernel %in% c("linear", "kernel"), "degree"] <- NA
      grid[grid$kernel %in% c("linear", "kernel"), "coef0"] <- NA
    }
    grid <- unique(grid)
  }

  return(grid)
}

#' Manage \code{mod} Arguments
#'
#' This function organizes \code{mod} arguments passed to \code{pl} functions.
#'
#' @param func A character string. The \code{mod} function to call.
#' @param ... Additional arguments passed to the \code{mod} function.
#' @return A list of arguments.
#' @export
ctrlModSet <- function(func, ...){

  list("func" = func, ...)
}

#' Manage \code{split} Arguments
#'
#' This function organizes \code{split} arguments passed to \code{pl} functions.
#'
#' @param func A character string. The \code{split} function to call.
#' @param percent.include Argument passed to the \code{split} function.
#' @param ... Additional arguments passed to the \code{split} function.
#' @return A list of arguments.
#' @export
ctrlSplitSet <- function(func, percent.include, ...){

  list("func" = func, "percent.include" = percent.include, ...)
}

#' Manage \code{fs} Arguments
#'
#' This function organizes \code{fs} arguments passed to \code{pl} functions.
#'
#' @param func A character string. The \code{fs} function to call.
#' @param top Argument passed to the \code{fs} function.
#' @param ... Additional arguments passed to the \code{fs} function.
#' @return A list of arguments.
#' @export
ctrlFeatureSelect <- function(func, top, ...){

  list("func" = func, "top" = top, ...)
}

#' Manage \code{plGrid} Arguments
#'
#' This function organizes \code{plGrid} arguments passed to \code{pl} functions.
#'
#' @param func A character string. The \code{pl} function to call.
#' @param top Argument passed to the \code{pl} function. Leave missing
#'  when handling \code{plMonteCarlo} or \code{plNested} arguments.
#' @param ... Additional arguments passed to the \code{pl} function.
#' @return A list of arguments.
#' @export
ctrlGridSearch <- function(func, top, ...){

  if(missing(top)){ list("func" = func, ...)
  }else{ list("func" = func, "top" = top, ...) }
}

#' Check \code{ctrlGS} Arguments
#'
#' This function ensures that the list of arguments for \code{ctrlGS} meets
#'  the criteria required by the \code{\link{plNested}} function. This
#'  function forces \code{aucSkip = TRUE} and \code{plotSkip = TRUE}.
#'
#' @param args A list of arguments to check.
check.ctrlGS <- function(args){

  if(args$func == "plGrid" | args$func == "plGridMulti"){
    if(!"aucSkip" %in% names(args)) args <- append(args, list("aucSkip" = TRUE))
    if(!args$aucSkip) args$aucSkip <- TRUE
    if(!"plotSkip" %in% names(args)) args <- append(args, list("plotSkip" = TRUE))
    if(!args$plotSkip) args$plotSkip <- TRUE
  }

  return(args)
}
