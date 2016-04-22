###########################################################
### Nested cross-validation argument check

#' Check \code{ctrlGS} Arguments
#'
#' This function ensures that a list of arguments for \code{ctrlGS} meets the
#'  criteria for the \code{\link{plNested}} function.
#'
#' This function mandates a non-NULL \code{fold} argument and a \code{TRUE}
#'  \code{aucSkip} argument.
#'
#' @param args A list of arguments to check.
check.ctrlGS <- function(args){

  if(args$func == "plGrid" & !"fold" %in% names(args)){

    cat("Setting 'fold' to 10 (default behavior, override explicitly)...\n")
    args <- append(args, list("fold" = 10))
  }

  if(args$func == "plGrid" & is.null(args$fold)){

    cat("Uh oh! This function requires non-NULL 'fold'. Setting 'fold' to 10...\n")
    args$fold <- 10
  }

  if(args$func == "plGrid" & !"aucSkip" %in% names(args)){

    cat("Setting 'aucSkip' to TRUE (default behavior, override explicitly)...\n")
    args <- append(args, list("aucSkip" = TRUE))
  }

  if(args$func == "plGrid" & !args$aucSkip){

    cat("Uh oh! This function requires TRUE 'aucSkip'. Setting 'aucSkip' to TRUE...\n")
    args$aucSkip <- TRUE
  }

  return(args)
}

###########################################################
### Nested cross-validation

#' Nested Cross-Validation
#'
#' Perform nested cross-validation.
#'
#' Analogous to how \code{\link{plGrid}} manages multiple \code{build} and
#'  \code{predict} tasks, \code{plNested} effectively manages multiple
#'  \code{plGrid} tasks.
#'
#' Specifically, \code{plNested} segregates the data into v-folds,
#'  treating each fold as a validation set to a training set built from
#'  the remaining subjects. Then, some \code{fold} times, it performs all
#'  feature selection tasks (listed via \code{ctrlFS}) on each split
#'  training set, and executing \code{plGrid} (via \code{ctrlGS}) on the
#'  training set.
#'
#' To perform multiple feature selection tasks, supply a list of multiple
#'  \code{\link{ctrlFeatureSelect}} argument wrappers to \code{ctrlFS}.
#'  To reduce the results of \code{plNested} to a single performance metric,
#'  feed the returned \code{ExprsPipeline} object through
#'  \code{\link{calcNested}}.
#'
#' Note that \code{plGrid} uses \code{\link{plCV}} to calculate the "inner"
#'  cross-validation accuracy. Depending on the use case, may not represent
#'  the most appropriate choice. We hope future implementations will accomodate
#'  other \code{ctrlGS} functions beyond \code{plGrid}.
#'
#' @param array Specifies the \code{ExprsArray} object to undergo cross-validation.
#' @param fold A numeric scalar. Specifies the number of folds for cross-validation.
#'  Set \code{fold = 0} to perform leave-one-out cross-validation.
#' @param ctrlFS A list of arguments handled by \code{\link{ctrlFeatureSelect}}.
#' @param ctrlGS Arguments handled by \code{\link{ctrlGridSearch}}.
#' @param save A logical scalar. Toggles whether to save each fold.
#' @return An \code{\link{ExprsPredict-class}} object.
#'
#' @seealso
#' \code{\link{plCV}}, \code{\link{plGrid}}, \code{\link{plMonteCarlo}}, \code{\link{plNested}}
#'
#' @examples
#' \dontrun{
#' require(golubEsets)
#' data(Golub_Merge)
#' array <- arrayEset(Golub_Merge, colBy = "ALL.AML", include = list("ALL", "AML"))
#' array <- modFilter(array, 20, 16000, 500, 5) # pre-filter Golub ala Deb 2003
#' array <- modTransform(array) # lg transform
#' array <- modNormalize(array, c(1, 2)) # normalize gene and subject vectors
#' fs <- ctrlFeatureSelect(func = "fsEbayes", probes = 0)
#' gs <- ctrlGridSearch(func = "plGrid", how = "buildANN", probes = c(10, 20, 30),
#'                      size = 1:3, decay = c(0, .5, 1), fold = 0)
#' nest <- plNested(arrays[[1]], fold = 10, ctrlFS = fs, ctrlGS = gs, save = FALSE)
#' }
#' @export
plNested <- function(array, fold = 10, ctrlFS, ctrlGS, save = FALSE){

  if(ctrlGS$func != "plGrid"){

    stop("Uh oh! This function currently only supports 'ctrlGS$func = plGrid'!")
  }

  if(!is.null(array@preFilter) | !is.null(array@reductionModel)){

    warning("Prior use of feature selection may result in ",
            "overly optimistic cross-validation accuracies!")
  }

  # Perform LOOCV if 0 fold
  if(fold == 0) fold <- nrow(array@annot)

  if(fold > nrow(array@annot)){

    warning("Insufficient subjects for v-fold cross-validation. Performing LOOCV instead.")
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

    # Perform fs_ function for each argument set in ctrlFS
    if(!"list" %in% lapply(ctrlFS, class)) ctrlFS <- list(ctrlFS)
    for(i in 1:length(ctrlFS)){

      func <- ctrlFS[[i]]$func
      args <- append(list("object" = array.boot), ctrlFS[[i]][!ctrlFS[[i]] %in% func])
      array.boot <- do.call(what = func, args = args)
    }

    # Perform some gridsearch function (e.g. plGrid)
    func <- ctrlGS$func
    args <- append(list("array.train" = array.boot,
                        "array.valid" = array.demi),
                   ctrlGS[!ctrlGS %in% func])
    args <- exprso:::checkArgs.ctrlGS(args)

    # Save pl object
    pl <- do.call(what = func, args = args)
    pl@summary <- cbind(v, pl@summary)
    pls[[v]] <- pl
  }

  pl <- new("ExprsPipeline",
            summary = do.call(rbind, lapply(pls, function(obj) obj@summary)),
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
#'  summarizing classifier performance. Listing multiple columns will calculate
#'  performance as a product of those listed performances.
#' @return A numeric scalar. The cross-validation accuracy.
#'
#' @export
calcNested <- function(pl, colBy){

  if(missing(colBy)) stop("Uh oh! Missing 'colBy' argument.")

  if("v" %in% colnames(pl@summary)){

    if("train.plCV" %in% colnames(pl@summary)){

      # Prepare container to store validation accuracy
      acc <- vector("numeric", length(unique(pl@summary$v)))

      for(b in 1:length(unique(pl@summary$v))){

        cat("Retrieving best accuracy for fold", b, "...\n")

        # Subset only fold 'b'
        fold <- pl@summary[pl@summary$v == b, ]

        # Select best model based on cross-validation accuracy
        best <- fold[which.max(fold$train.plCV), ]

        # Save validation accuracy as colBy product
        acc[b] <- apply(best[colBy], MARGIN = 1, prod)
      }

    }else{

      stop("Uh oh! Supplied data not in expected format. ",
           "Cannot calculate this cross-validation accuracy.")
    }

  }else{

    stop("Uh oh! Supplied data not in expected format. ",
         "Cannot calculate this cross-validation accuracy.")
  }

  # Return average validation accuracy
  cat("Averaging best accuracies across all folds...\n")
  return(mean(acc))
}
