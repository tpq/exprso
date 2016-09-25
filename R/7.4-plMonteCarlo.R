###########################################################
### Cross-validation argument handlers

#' Manage \code{split} Arguments
#'
#' This function organizes \code{split} arguments passed to \code{pl} functions.
#'
#' @param func A character string. The \code{split} function to call.
#' @param percent.include Argument passed to the \code{split} function.
#' @param ... Additional arguments passed to the \code{split} function.
#' @return A list of arguments.
#'
#' @export
ctrlSplitSet <- function(func, percent.include, ...){

  list("func" = func,
       "percent.include" = percent.include,
       ...)
}

#' Manage \code{fs} Arguments
#'
#' This function organizes \code{fs} arguments passed to \code{pl} functions.
#'
#' @param func A character string. The \code{fs} function to call.
#' @param top Argument passed to the \code{fs} function.
#' @param ... Additional arguments passed to the \code{fs} function.
#' @return A list of arguments.
#'
#' @export
ctrlFeatureSelect <- function(func, top, ...){

  list("func" = func,
       "top" = top,
       ...)
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
#'
#' @export
ctrlGridSearch <- function(func, top, ...){

  if(missing(top)){

    list("func" = func,
         ...)
  }else{

    list("func" = func,
         "top" = top,
         ...)
  }
}

###########################################################
### Monte Carlo cross-validation

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
#' @param array Specifies the \code{ExprsArray} object to undergo cross-validation.
#' @param B A numeric scalar. The number of times to \code{split} the data.
#' @param ctrlSS Arguments handled by \code{\link{ctrlSplitSet}}.
#' @param ctrlFS A list of arguments handled by \code{\link{ctrlFeatureSelect}}.
#' @param ctrlGS Arguments handled by \code{\link{ctrlGridSearch}}.
#' @param save A logical scalar. Toggles whether to save randomly split
#'  training and validation sets.
#'
#' @return An \code{\link{ExprsPipeline-class}} object.
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
plMonteCarlo <- function(array, B = 10, ctrlSS, ctrlFS, ctrlGS, save = FALSE){

  # For each bootstrap
  pls <- lapply(1:B,
                function(boot){

                  # Perform some split function (e.g. splitStrat)
                  func <- ctrlSS$func
                  args <- append(list("object" = array), ctrlSS[!ctrlSS %in% func])
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

                  # Perform some gridsearch function (e.g. plGrid)
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

                  }else{

                    stop("Uh oh! No method in place for this 'pl' pipeline.")
                  }


                  # Append pl@summary
                  pl@summary <- cbind(boot, pl@summary)

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
#'  summarizing classifier performance. Listing multiple columns will calculate
#'  performance as a product of those listed performances.
#' @return A numeric scalar. The cross-validation accuracy.
#'
#' @export
calcMonteCarlo <- function(pl, colBy){

  if(missing(colBy)) stop("Uh oh! Missing 'colBy' argument.")

  if("boot" %in% colnames(pl@summary)){

    if("train.plCV" %in% colnames(pl@summary)){

      # Prepare container to store validation accuracy
      acc <- vector("numeric", length(unique(pl@summary$boot)))

      for(b in 1:length(unique(pl@summary$boot))){

        cat("Retrieving best accuracy for boot", b, "...\n")

        # Subset only boot 'b'
        boot <- pl@summary[pl@summary$boot == b, ]

        # Select best model based on cross-validation accuracy
        best <- boot[which.max(boot$train.plCV), ]

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
  cat("Averaging best accuracies across all boots...\n")
  return(mean(acc))
}
