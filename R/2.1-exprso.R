#' The \code{exprso} Package
#'
#' @description
#' Welcome to the \code{exprso} package!
#'
#' The \code{exprso} function imports data into the learning environment.
#'
#' See \code{\link{mod}} to pre-process the data.
#'
#' See \code{\link{split}} to split off a test set.
#'
#' See \code{\link{fs}} to select features.
#'
#' See \code{\link{build}} to build models.
#'
#' See \code{\link{pl}} to build models high-throughput.
#'
#' See \code{\link{buildEnsemble}} to build ensembles.
#'
#' See \code{\link{pipeFilter}} to pre-process ensembles.
#'
#' See \code{\link{exprso-predict}} to deply models.
#'
#' @param x A matrix of feature data for all samples. Rows should
#'  contain samples and columns should contain features.
#' @param y A vector of outcomes for all samples. If
#'  \code{class(y) == "character"} or \code{class(y) == "factor"},
#'  \code{exprso} prepares data for binary or multi-class classification.
#'  Else, \code{exprso} prepares data for regression. If \code{y} is a
#'  matrix, the program assumes the first column is the outcome.
#' @return An \code{ExprsArray} object.
#' @export
exprso <- function(x, y){

  if(length(y) != nrow(x)) stop("Incorrect number of outcomes.")
  array <-
    new("ExprsArray",
        exprs = t(as.data.frame(x)), annot = as.data.frame(y),
        preFilter = NULL, reductionModel = NULL
    )

  # Prepare ExprsArray object using x and y input
  colnames(array@exprs) <- paste0("x", 1:ncol(array@exprs))
  colnames(array@exprs) <- make.names(colnames(array@exprs), unique = TRUE)
  rownames(array@annot) <- colnames(array@exprs)
  labels <- array@annot[,1]

  # Set sub-class to guide fs and build modules
  if(class(labels) == "character" | class(labels) == "factor"){
    if(length(unique(y)) == 2){
      print("Preparing data for binary classification.")
      class(array) <- "ExprsBinary"
      array@annot$defineCase <- ifelse(labels == unique(labels)[1], "Control", "Case")
    }else{
      print("Preparing data for multi-class classification.")
      class(array) <- "ExprsMulti"
      array@annot$defineCase <- factor(labels)
    }
  }else{
    print("Preparing data for regression.")
    class(array) <- "ExprsCont"
    array@annot$defineCase <- labels
  }

  # Remove features with any NA values
  if(any(is.na(array@exprs))){
    print("Removing features with NA values.")
    noNAs <- apply(array@exprs, 1, function(x) !any(is.na(x)))
    array@exprs <- array@exprs[noNAs, ]
  }

  return(array)
}

#' @name mod
#' @rdname mod
#'
#' @title Process Data
#'
#' @description
#' The \code{exprso} package includes these pre-process modules:
#'
#' - \code{\link{modFilter}}
#'
#' - \code{\link{modTransform}}
#'
#' - \code{\link{modNormalize}}
#'
#' - \code{\link{modTMM}}
NULL

#' @name pl
#' @rdname pl
#'
#' @title Classification Pipelines
#'
#' @description
#' The \code{exprso} package includes these automated pipeline modules:
#'
#' - \code{\link{plCV}}
#'
#' - \code{\link{plGrid}}
#'
#' - \code{\link{plGridMulti}}
#'
#' - \code{\link{plMonteCarlo}}
#'
#' - \code{\link{plNested}}
NULL
