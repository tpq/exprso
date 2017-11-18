###########################################################
### Simple cross-validation

#' Perform Simple Cross-Validation
#'
#' Calculates v-fold or leave-one-out cross-validation without selecting a new
#'  set of features with each fold.
#'
#' \code{plCV} performs v-fold or leave-one-out cross-validation for an
#'  \code{ExprsArray} object. The argument \code{fold} specifies the number
#'  of v-folds to use during cross-validation. Set \code{fold = 0} to perform
#'  leave-one-out cross-validation. This approach to cross-validation
#'  will work for \code{ExprsBinary} and \code{ExprsMulti} objects alike. The
#'  performance metric used to measure cross-validation accuracy is the
#'  \code{acc} slot returned by \code{\link{calcStats}}.
#'
#' This type of cross-validation is most appropriate if the \code{ExprsArray}
#'  has not undergone any prior feature selection. However, it may also have a role
#'  as an unbiased guide to parameter selection when embedded in
#'  \code{\link{plGrid}}. If using cross-validation in lieu of an independent test
#'  set in the setting of one or more feature selection methods, consider using
#'  a more "sophisticated" form of cross-validation as implemented in
#'  \code{\link{plMonteCarlo}} or \code{\link{plNested}}.
#'
#' When calculating classifier performance with \code{\link{calcStats}}, this
#'  function forces \code{aucSkip = TRUE} and \code{plotSkip = TRUE}.
#'
#' @param array Specifies the \code{ExprsArray} object to undergo cross-validation.
#' @inheritParams fs.
#' @param how A character string. Specifies the \code{\link{build}} method to iterate.
#' @param fold A numeric scalar. Specifies the number of folds for cross-validation.
#'  Set \code{fold = 0} to perform leave-one-out cross-validation.
#' @param ... Arguments passed to the \code{how} method.
#'
#' @return A numeric scalar. The cross-validation accuracy.
#'
#' @seealso
#' \code{\link{fs}}\cr
#' \code{\link{build}}\cr
#' \code{\link{doMulti}}\cr
#' \code{\link{exprso-predict}}\cr
#' \code{\link{plCV}}\cr
#' \code{\link{plGrid}}\cr
#' \code{\link{plGridMulti}}\cr
#' \code{\link{plMonteCarlo}}\cr
#' \code{\link{plNested}}
#'
#' @export
plCV <- function(array, top, how, fold, ...){

  args <- as.list(substitute(list(...)))[-1]

  # Perform LOOCV if 0 fold
  if(fold == 0) fold <- nrow(array@annot)

  if(fold > nrow(array@annot)){

    warning("Insufficient subjects for plCV v-fold cross-validation. Performing LOOCV instead.")
    fold <- nrow(array@annot)
  }

  # Add the ith subject ID to the vth fold
  subjects <- vector("list", fold)
  ids <- sample(rownames(array@annot))
  i <- 1
  while(i <= nrow(array@annot)){

    subjects[[i %% fold + 1]] <- c(subjects[[i %% fold + 1]], ids[i])
    i <- i + 1
  }

  # Build a machine against the vth fold
  accs <- vector("numeric", fold)
  for(v in 1:length(subjects)){

    # The leave one out
    array.train <- new(class(array),
                       exprs = array@exprs[, !colnames(array@exprs) %in% subjects[[v]], drop = FALSE],
                       annot = array@annot[!rownames(array@annot) %in% subjects[[v]], ],
                       preFilter = array@preFilter,
                       reductionModel = array@reductionModel)

    # The left out one
    array.valid <- new(class(array),
                       exprs = array@exprs[, colnames(array@exprs) %in% subjects[[v]], drop = FALSE],
                       annot = array@annot[rownames(array@annot) %in% subjects[[v]], ],
                       preFilter = array@preFilter,
                       reductionModel = array@reductionModel)

    # Prepare args for do.call
    args.v <- append(list("object" = array.train, "top" = top), args)

    # Build machine and deploy
    mach <- do.call(what = how, args = args.v)
    pred <- predict(mach, array.valid, verbose = FALSE)
    accs[v] <- calcStats(pred, aucSkip = TRUE, plotSkip = TRUE)$acc

    cat("plCV", v, "accuracy:", accs[v], "\n")
  }

  acc <- mean(accs)

  return(acc)
}
