#' Nested Cross-Validation
#'
#' Perform nested cross-validation.
#'
#' Analogous to how \code{\link{plGrid}} manages multiple \code{build} and
#'  \code{predict} tasks, one can think of \code{plNested} as managing
#'  multiple \code{pl} tasks.
#'
#' Specifically, \code{plNested} segregates the data into v-folds,
#'  treating each fold as a validation set and the subjects not in that fold
#'  as a training set. Then, some \code{fold} times, it performs all
#'  feature selection tasks (listed via \code{ctrlFS}) on each split
#'  of the data, and executes the \code{pl} function (via \code{ctrlGS})
#'  using the training set.
#'
#' To perform multiple feature selection tasks, supply a list of multiple
#'  \code{\link{ctrlFeatureSelect}} argument wrappers to \code{ctrlFS}.
#'  To reduce the results of \code{plNested} to a single performance metric,
#'  you can feed the returned \code{ExprsPipeline} object through the helper
#'  function \code{\link{calcNested}}.
#'
#' When calculating model performance with \code{\link{calcStats}}, this
#'  function forces \code{aucSkip = TRUE} and \code{plotSkip = TRUE}.
#'  When embedding another \code{plMonteCarlo} or \code{plNested} call within
#'  this function (i.e., via \code{ctrlGS}), outer-fold model performance
#'  will force \code{aucSkip = TRUE} and \code{plotSkip = TRUE}.
#'
#' @param array Specifies the \code{ExprsArray} object to undergo cross-validation.
#' @param fold A numeric scalar. Specifies the number of folds for cross-validation.
#'  Set \code{fold = 0} to perform leave-one-out cross-validation.
#' @param ctrlFS A list of arguments handled by \code{\link{ctrlFeatureSelect}}.
#' @param ctrlGS Arguments handled by \code{\link{ctrlGridSearch}}.
#' @param save A logical scalar. Toggles whether to save each fold.
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
#' fs <- ctrlFeatureSelect(func = "fsEbayes", top = 0)
#' gs <- ctrlGridSearch(func = "plGrid", how = "buildANN", top = c(10, 20, 30),
#'                      size = 1:3, decay = c(0, .5, 1), fold = 0)
#' nest <- plNested(arrays[[1]], fold = 10, ctrlFS = fs, ctrlGS = gs, save = FALSE)
#' }
#' @export
plNested <- function(array, fold = 10, ctrlFS = NULL, ctrlGS, save = FALSE){

  # Perform LOOCV if 0 fold
  if(fold == 0) fold <- nrow(array@annot)

  if(fold > nrow(array@annot)){

    warning("Insufficient subjects for plNested v-fold cross-validation. Performing LOOCV instead.")
    fold <- nrow(array@annot)
  }

  # Add the ith random subject ID to the vth fold
  subjects <- vector("list", fold)
  ids <- sample(rownames(array@annot))
  i <- 1
  while(i <= nrow(array@annot)){

    subjects[[i %% fold + 1]] <- c(subjects[[i %% fold + 1]], ids[i])
    i <- i + 1
  }

  # Perform nested cross-validation
  pls <- vector("list", fold)
  numTicks <- 0
  for(v in 1:length(subjects)){

    # The v-th fold
    array.boot <- new(class(array),
                      exprs = array@exprs[, !colnames(array@exprs) %in% subjects[[v]], drop = FALSE],
                      annot = array@annot[!rownames(array@annot) %in% subjects[[v]], ],
                      preFilter = array@preFilter,
                      reductionModel = array@reductionModel)

    # The leave out
    array.demi <- new(class(array),
                      exprs = array@exprs[, colnames(array@exprs) %in% subjects[[v]], drop = FALSE],
                      annot = array@annot[rownames(array@annot) %in% subjects[[v]], ],
                      preFilter = array@preFilter,
                      reductionModel = array@reductionModel)

    # Save files
    if(save){

      save(array.boot, file = paste0("plNested ", v, " (",
                                     gsub(":", ".", Sys.time()),
                                     ") bootstrap.RData"))
      save(array.demi, file = paste0("plNested ", v, " (",
                                     gsub(":", ".", Sys.time()),
                                     ") demi-holdout.RData"))
    }

    # Optionally perform fs_ function for each argument set in ctrlFS
    if(!is.null(ctrlFS)){

      if(!"list" %in% lapply(ctrlFS, class)) ctrlFS <- list(ctrlFS)
      for(i in 1:length(ctrlFS)){

        func <- ctrlFS[[i]]$func
        args <- append(list("object" = array.boot), ctrlFS[[i]][!ctrlFS[[i]] %in% func])
        array.boot <- suppressMessages(do.call(what = func, args = args))
      }
    }

    # Perform some gridsearch function (e.g., plGrid)
    if(ctrlGS$func %in% c("plGrid", "plGridMulti")){

      # Update progress bar for inner plGrid gridsearch
      numTicks <- progress(v, fold, numTicks)

      args <- append(list("array.train" = array.boot,
                          "array.valid" = array.demi),
                     ctrlGS)
      args <- check.ctrlGS(args)
      func <- ctrlGS$func
      pl <- suppressMessages(do.call(what = func, args = args[!args %in% func]))

    }else if(ctrlGS$func %in% c("plMonteCarlo", "plNested")){

      args <- append(list("array" = array.boot), ctrlGS)
      func <- ctrlGS$func
      pl <- do.call(what = func, args = args[!args %in% func])

      # Calculate outer fold accuracy for embedded pipelines
      preds <- sapply(pl@machs, predict, array.demi)
      stats <- lapply(preds, calcStats, aucSkip = TRUE, plotSkip = TRUE)
      pl@summary <- cbind("outer" = do.call(rbind, stats), pl@summary)

    }else{

      stop("Uh oh! No method in place for this 'pl' pipeline.")
    }

    # Append pl@summary
    pl@summary <- cbind("pl" = "plNested", v,
                        "ss" = paste0(fold, "-fold"),
                        "fs" = paste(sapply(ctrlFS, "[", "func"),
                                     collapse = ", "),
                        "gs" = ctrlGS$func,
                        pl@summary)

    colnames(pl@summary) <- make.names(colnames(pl@summary), unique = TRUE)
    pls[[v]] <- pl
  }

  pl <- new("ExprsPipeline",
            summary = do.call(plyr::rbind.fill, lapply(pls, function(obj) obj@summary)),
            machs = unlist(lapply(pls, function(obj) obj@machs))
  )
}

#' Calculate \code{plNested} Performance
#'
#' \code{calcNested} calculates a single performance measure for the
#'  results of a \code{plNested} function call.
#'
#' For each dataset split (i.e., bootstrap), \code{calcNested} averages
#'  the validation set performance for the "best" model (where "best" is
#'  defined as the model with the maximum "internal" cross-validation
#'  accuracy, \code{max($train.plCV)}). The validation set performance
#'  ultimately averaged depends on the supplied \code{colBy} argument.
#'
#' @param pl Specifies the \code{ExprsPipeline} object returned by \code{plNested}.
#' @param colBy A character vector or string. Specifies column(s) to use when
#'  summarizing model performance. Listing multiple columns will calculate
#'  performance as a product of those listed performances.
#' @return A numeric scalar. The cross-validation accuracy.
#'
#' @export
calcNested <- function(pl, colBy = "valid.acc"){

  if(!"v" %in% colnames(pl@summary) | !"train.plCV" %in% colnames(pl@summary)){
    stop("Summary must have 'v' and 'train.plCV' columns.")
  }

  acc <- vector("numeric", length(unique(pl@summary$v)))
  for(b in 1:length(unique(pl@summary$v))){

    # Subset only fold 'b'
    fold <- pl@summary[pl@summary$v == b, ]

    # Select best model based on cross-validation accuracy
    best <- fold[which.max(fold$train.plCV), ]

    # Save validation accuracy as colBy product
    acc[b] <- apply(best[colBy], MARGIN = 1, prod)
  }

  # Return average validation accuracy
  cat("Averaging best accuracies across all folds...\n")
  return(mean(acc))
}
