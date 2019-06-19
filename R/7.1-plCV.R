#' Perform Simple Cross-Validation
#'
#' Calculates v-fold or leave-one-out cross-validation without selecting a new
#'  set of features with each fold. See Details.
#'
#' \code{plCV} performs v-fold or leave-one-out cross-validation. The argument
#'  \code{fold} specifies the number of v-folds to use during cross-validation.
#'  Set \code{fold = 0} to perform leave-one-out cross-validation.
#'
#' This type of cross-validation is most appropriate if the data
#'  has not undergone any prior feature selection. However, it is also useful
#'  as an unbiased guide to parameter selection within another
#'  \code{\link{pl}} workflow.
#'
#' Users should never need to call this function directly. Instead, they
#'  should use \code{\link{plMonteCarlo}} or \code{\link{plNested}}.
#'  There, \code{plCV} handles inner-fold cross-validation.
#'
#' @param array Specifies the \code{ExprsArray} object to undergo cross-validation.
#' @inheritParams plGrid
#' @return The average inner-fold cross-validation accuracy.
#' @export
plCV <- function(array, top, how, fold, aucSkip, plCV.acc, ...){

  args.how <- getArgs(...)

  # Perform LOOCV if 0 fold
  if(fold == 0) fold <- nsamps(array)
  if(fold > nsamps(array)){

    warning("Insufficient subjects for plCV v-fold cross-validation. Performing LOOCV instead.")
    fold <- nsamps(array)
  }

  # Add the ith subject ID to the vth fold
  ids <- sample(rownames(array@annot))
  splits <- suppressWarnings(split(ids, 1:fold)) # warns that some splits are bigger than others

  # Build a machine against the vth fold
  accs <- vector("numeric", fold)
  for(v in 1:length(splits)){

    holdout <- colnames(array@exprs) %in% splits[[v]]
    array.train <- array[!holdout,]
    array.valid <- array[holdout,]

    # Build machine and deploy
    args.v <- append(list("object" = array.train, "top" = top), args.how)
    mach.v <- do.call(what = how, args = args.v)
    pred.v <- predict(mach.v, array.valid, verbose = FALSE)
    accs[v] <- calcStats(pred.v, aucSkip = aucSkip, plotSkip = TRUE, verbose = FALSE)[, plCV.acc]
  }

  acc <- mean(accs)

  return(acc)
}
