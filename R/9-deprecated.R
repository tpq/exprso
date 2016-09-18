###########################################################
### Import data from file or eSet

#' Import Data from File
#'
#' A convenience function that builds an \code{ExprsArray} object from a tab-delimited file.
#'
#' \code{arrayRead} helps build an \code{ExprsArray} object from a tab-delimited
#'  data file, passing along the \code{file} and \code{...} argument(s) to
#'  \code{\link{read.delim}}. This function expects that the delimited data file has
#'  the following format: rows indicate subject entries while columns indicate measured
#'  variables.
#'
#' The first several columns should contain annotation information (e.g., age,
#'  sex, diagnosis). The remaining columns should contain feature data (e.g. expression
#'  values). The argument \code{probes.begin} defines the j-th column at which the feature
#'  data starts. By default, \code{arrayRead} forces \code{stringsAsFactors = FASE}.
#'
#' This function automatically removes any features with \code{NA} values.
#'
#' @param file A character string. Argument passed along to \code{read.delim}.
#' @param probes.begin A numeric scalar. The j-th column at which feature data starts.
#' @param colID A numeric or character index. The column used to name subjects.
#' @param colBy A numeric or character index. The column that contains group annotations.
#' @param include A list of character vectors. Specifies which annotations in \code{colBy}
#'  to include into which groups. Each element of a list specifies a unique group while
#'  each element of the character vector specifies an annotation to fit to that group. For
#'  binary classification, the first list element defines the negative or control group.
#' @param ... Additional arguments passed along to \code{read.delim}.
#'
#' @seealso
#' \code{\link{ExprsArray-class}}, \code{\link{arrayEset}}, \code{\link{GSE2eSet}}
#' @importFrom utils read.delim
#' @export
arrayRead <- function(file, probes.begin, colID, colBy, include, ...){

  warning("This function is now deprecated. Use 'arrayExprs' from here-on!")

  if(!class(include) == "list") stop("Uh oh! User must provide 'include' argument as list!")
  if(length(include) < 2) stop("Uh oh! User must provide at least two classes!")

  # Pass any additional arguments to read.delim
  table <- read.delim(file, stringsAsFactors = FALSE, ...)

  # Label table by proper subject ID names
  rownames(table) <- make.names(table[, colID], unique = TRUE)

  # Separate the expression data from annotation data
  array <- new(ifelse(length(include) == 2, "ExprsBinary", "ExprsMulti"),
               exprs = t(table[, probes.begin:ncol(table)]),
               annot = table[, 1:(probes.begin-1)],
               preFilter = NULL,
               reductionModel = NULL)

  # Filter out samples in @annot not in 'include'
  array@annot <- array@annot[array@annot[, colBy] %in% unlist(include), ]

  # Number classes in order of 'include'
  for(i in 1:length(include)){

    array@annot[array@annot[, colBy] %in% include[[i]], "defineCase"] <- i
  }

  # Name classes if binary
  if(class(array) == "ExprsBinary"){

    array@annot$defineCase <- ifelse(array@annot$defineCase == 1, "Control", "Case")
  }

  # Set factor if multi
  if(class(array) == "ExprsMulti"){

    array@annot$defineCase <- factor(array@annot$defineCase)
  }

  # Filter @exprs
  array@exprs <- array@exprs[, rownames(array@annot)]

  # Remove probes with missing values
  if(any(is.na(array@exprs))){

    cat("Removing probes with missing values...\n")
    array@exprs <- array@exprs[apply(array@exprs, 1, function(x) !any(is.na(x))), ]
  }

  return(array)
}

#' Import Data from eSet
#'
#' A convenience function that builds an \code{ExprsArray} object from an \code{eSet} object.
#'
#' The package Biobase maintains a popular class object called \code{ExpressionSet} that
#'  often gets used to store expression data. This function converts this \code{eSet}
#'  object into an \code{ExprsArray} object. This function automatically removes any
#'  features with \code{NA} values.
#'
#' @param eSet An \code{eSet} object.
#' @inheritParams arrayRead
#'
#' @seealso
#' \code{\link{ExprsArray-class}}, \code{\link{arrayRead}}, \code{\link{GSE2eSet}}
#' @import Biobase
#' @export
arrayEset <- function(eSet, colBy, include){

  if(!class(include) == "list") stop("Uh oh! User must provide 'include' argument as list!")
  if(length(include) < 2) stop("Uh oh! User must provide at least two classes!")

  # Build an ExprsArray object from the eSet object
  array <- new(ifelse(length(include) == 2, "ExprsBinary", "ExprsMulti"),
               exprs = exprs(eSet),
               annot = eSet@phenoData@data,
               preFilter = NULL,
               reductionModel = NULL)

  # Force @annot rownames to mirror proper @exprs colnames
  colnames(array@exprs) <- make.names(colnames(array@exprs), unique = TRUE)
  rownames(array@annot) <- colnames(array@exprs)

  # Filter out samples in @annot not in 'include'
  array@annot <- array@annot[array@annot[, colBy] %in% unlist(include), ]

  # Number classes in order of 'include'
  for(i in 1:length(include)){

    array@annot[array@annot[, colBy] %in% include[[i]], "defineCase"] <- i
  }

  # Name classes if binary
  if(class(array) == "ExprsBinary"){

    array@annot$defineCase <- ifelse(array@annot$defineCase == 1, "Control", "Case")
  }

  # Set factor if multi
  if(class(array) == "ExprsMulti"){

    array@annot$defineCase <- factor(array@annot$defineCase)
  }

  # Filter @exprs
  array@exprs <- array@exprs[, rownames(array@annot)]

  # Remove probes with missing values
  if(any(is.na(array@exprs))){

    cat("Removing probes with missing values...\n")
    array@exprs <- array@exprs[apply(array@exprs, 1, function(x) !any(is.na(x))), ]
  }

  return(array)
}
