#' Monte Carlo Cross-Validation
#'
#' Perform Monte Carlo style cross-validation.
#'
#' Analogous to how \code{\link{plGrid}} manages multiple \code{build} and
#'  \code{predict} tasks, one can think of \code{plMonteCarlo} as managing
#'  multiple \code{pl} tasks.
#'
#' Specifically, \code{plMonteCarlo} will call the provided \code{split}
#'  function (via \code{ctrlSS}) some \code{B} times, perform all
#'  feature selection tasks (listed via \code{ctrlFS}) on each split of
#'  the data, and execute the \code{pl} function (via \code{ctrlGS})
#'  using the bootstrapped set.
#'
#' To perform multiple feature selection tasks, supply a list of multiple
#'  \code{\link{ctrlFeatureSelect}} argument wrappers to \code{ctrlFS}.
#'  To reduce the results of \code{plMonteCarlo} to a single performance metric,
#'  you can feed the returned \code{ExprsPipeline} object through the helper
#'  function \code{\link{calcMonteCarlo}}.
#'
#' When embedding another \code{plMonteCarlo} or \code{plNested} call within
#'  this function (i.e., via \code{ctrlGS}), outer-fold model performance
#'  will force \code{aucSkip = TRUE} and \code{plotSkip = TRUE}.
#'
#' @param array Specifies the \code{ExprsArray} object to undergo cross-validation.
#' @param B A numeric scalar. The number of times to \code{split} the data.
#' @param ctrlSS Arguments handled by \code{\link{ctrlSplitSet}}.
#' @param ctrlFS A list of arguments handled by \code{\link{ctrlFeatureSelect}}.
#' @param ctrlGS Arguments handled by \code{\link{ctrlGridSearch}}.
#' @param ctrlMS Arguments handled by \code{\link{ctrlModSet}}. Optional.
#' @param save A logical scalar. Toggles whether to save randomly split
#'  training and validation sets.
#' @return An \code{\link{ExprsPipeline-class}} object.
#'
#' @examples
#' \dontrun{
#' require(golubEsets)
#' data(Golub_Merge)
#' array <- arrayEset(Golub_Merge, colBy = "ALL.AML", include = list("ALL", "AML"))
#' array <- modFilter(array, 20, 16000, 500, 5) # pre-filter Golub ala Deb 2003
#' array <- modTransform(array) # lg transform
#' array <- modNormalize(array, c(1, 2)) # normalize gene and subject vectors
#' ss <- ctrlSplitSet(func = "splitStratify", percent.include = 67, colBy = NULL)
#' fs <- list(ctrlFeatureSelect(func = "fsStats", top = 0, how = "t.test"),
#'            ctrlFeatureSelect(func = "fsPrcomp", top = 50))
#' gs <- ctrlGridSearch(func = "plGrid", how = "buildSVM", top = c(2, 3, 4), fold = 10,
#'                      kernel = c("linear", "radial"), cost = 10^(-3:3), gamma = 10^(-3:3))
#' boot <- plMonteCarlo(array, B = 3, ctrlSS = ss, ctrlFS = fs, ctrlGS = gs)
#' }
#' @export
plMonteCarlo <- function(array, B = 10, ctrlSS, ctrlFS, ctrlGS, ctrlMS = NULL, save = FALSE){

  # For each bootstrap
  pls <- lapply(1:B,
                function(boot){

                  # Optionally apply a mod function here
                  array.b <- array
                  if(!is.null(ctrlMS)){
                    func <- ctrlMS$func
                    args <- append(list("object" = array), ctrlMS[!ctrlMS %in% func])
                    array.b <- do.call(what = func, args = args)
                  }

                  # Perform some split function (e.g. splitStrat)
                  func <- ctrlSS$func
                  args <- append(list("object" = array.b), ctrlSS[!ctrlSS %in% func])
                  arrays <- do.call(what = func, args = args)
                  array.boot <- arrays[[1]]
                  array.demi <- arrays[[2]]

                  # Save files
                  if(save){

                    save(array.boot, file = paste0("plMonteCarlo ", boot, " (",
                                                   gsub(":", ".", Sys.time()),
                                                   ") bootstrap.RData"))

                    save(array.demi, file = paste0("plMonteCarlo ", boot, " (",
                                                   gsub(":", ".", Sys.time()),
                                                   ") demi-holdout.RData"))
                  }

                  # Perform fs_ function for each argument set in ctrlFS
                  if(!"list" %in% lapply(ctrlFS, class)) ctrlFS <- list(ctrlFS)
                  for(i in 1:length(ctrlFS)){

                    func <- ctrlFS[[i]]$func
                    args <- append(list("object" = array.boot), ctrlFS[[i]][!ctrlFS[[i]] %in% func])
                    array.boot <- do.call(what = func, args = args)
                  }

                  # Perform some gridsearch function (e.g., plGrid)
                  if(ctrlGS$func %in% c("plGrid", "plGridMulti")){

                    func <- ctrlGS$func
                    args <- append(list("array.train" = array.boot,
                                        "array.valid" = array.demi),
                                   ctrlGS[!ctrlGS %in% func])
                    pl <- do.call(what = func, args = args)

                  }else if(ctrlGS$func %in% c("plMonteCarlo", "plNested")){

                    func <- ctrlGS$func
                    args <- append(list("array" = array.boot), ctrlGS[!ctrlGS %in% func])
                    pl <- do.call(what = func, args = args)

                    # Calculate outer fold accuracy for embedded pipelines
                    preds <- sapply(pl@machs, predict, array.demi)
                    stats <- lapply(preds, calcStats, aucSkip = TRUE, plotSkip = TRUE)
                    pl@summary <- cbind("outer" = do.call(rbind, stats), pl@summary)

                  }else{

                    stop("Uh oh! No method in place for this 'pl' pipeline.")
                  }

                  # Append pl@summary
                  pl@summary <- cbind("pl" = "plMonteCarlo", boot,
                                      "ss" = ctrlSS$func,
                                      "fs" = paste(sapply(ctrlFS, "[", "func"),
                                                   collapse = ", "),
                                      "gs" = ctrlGS$func,
                                      pl@summary)

                  colnames(pl@summary) <- make.names(colnames(pl@summary), unique = TRUE)
                  return(pl)
                }
  )

  pl <- new("ExprsPipeline",
            summary = do.call(plyr::rbind.fill, lapply(pls, function(obj) obj@summary)),
            machs = unlist(lapply(pls, function(obj) obj@machs))
  )

  return(pl)
}

#' Calculate \code{plMonteCarlo} Performance
#'
#' \code{calcMonteCarlo} calculates a single performance measure for the
#'  results of a \code{plMonteCarlo} function call.
#'
#' For each dataset split (i.e., bootstrap), \code{calcMonteCarlo} averages
#'  the validation set performance for the "best" model (where "best" is
#'  defined as the model with the maximum "internal" cross-validation
#'  accuracy, \code{max($train.plCV)}). The validation set performance
#'  ultimately averaged depends on the supplied \code{colBy} argument.
#'
#' @param pl Specifies the \code{ExprsPipeline} object returned by \code{plMonteCarlo}.
#' @param colBy A character vector or string. Specifies column(s) to use when
#'  summarizing model performance. Listing multiple columns will calculate
#'  performance as a product of those listed performances.
#' @return A numeric scalar. The cross-validation accuracy.
#'
#' @export
calcMonteCarlo <- function(pl, colBy = "valid.acc"){

  if(!"boot" %in% colnames(pl@summary) | !"train.plCV" %in% colnames(pl@summary)){
    stop("Summary must have 'boot' and 'train.plCV' columns.")
  }

  acc <- vector("numeric", length(unique(pl@summary$boot)))
  for(b in 1:length(unique(pl@summary$boot))){

    # Subset only boot 'b'
    boot <- pl@summary[pl@summary$boot == b, ]

    # Select best model based on cross-validation accuracy
    best <- boot[which.max(boot$train.plCV), ]

    # Save validation accuracy as colBy product
    acc[b] <- apply(best[colBy], MARGIN = 1, prod)
  }

  # Return average validation accuracy
  cat("Averaging best accuracies across all boots...\n")
  return(mean(acc))
}
