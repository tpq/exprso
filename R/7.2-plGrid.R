###########################################################
### High-throughput classification

#' Perform High-Throughput Classification
#'
#' Trains and deploys multiple classifiers across a vast parameter search space.
#'
#' \code{plGrid} will \code{\link{build}} and \code{\link{exprso-predict}} for
#'  each combination of parameters listed as additional arguments (\code{...}).
#'  When using \code{plGrid}, supplying a numeric vector as the \code{probes}
#'  argument will train and deploy a classifier of each mentioned size for
#'  each combination of parameters listed in \code{...}. To skip prediction,
#'  set \code{array.valid = NULL}. Either way, this function returns an
#'  \code{\link{ExprsPipeline-class}} object which contains the build parameters,
#'  classifier performances, and \code{ExprsModel} objects for all trained
#'  and deployed models.
#'
#' \code{plGrid} will perform v-fold or leave-one-out cross-validation for the
#'  using \code{\link{plCV}}. The argument \code{fold} specifies the number
#'  of v-folds to use during cross-validation. Set \code{fold = 0} to perform
#'  leave-one-out cross-validation. This approach to cross-validation
#'  will work for \code{ExprsBinary} and \code{ExprsMulti} objects alike. The
#'  peformance metric used to measure cross-validation accuracy is the
#'  \code{acc} slot returned by \code{\link{calcStats}}. Set \code{fold = NULL}
#'  to skip cross-validation altogether.
#'
#' The use of \code{\link{plCV}} is most appropriate if the \code{ExprsArray}
#'  has not undergone any prior feature selection. However, it may also have a role
#'  as an unbiased guide to parameter selection when embedded in
#'  \code{\link{plGrid}}. If using cross-validation in lieu of an independent test
#'  set in the setting of one or more feature selection methods, consider using
#'  a more "sophisticated" form of cross-validation as implemented in
#'  \code{\link{plMonteCarlo}} or \code{\link{plNested}}.
#'
#' @param array.train Specifies the \code{ExprsArray} object to use as training set.
#' @param array.valid Specifies the \code{ExprsArray} object to use as validation set.
#' @param probes A numeric scalar or character vector. A numeric scalar indicates
#'  the number of top features that should undergo feature selection. A character vector
#'  indicates specifically which features by name should undergo feature selection.
#'  Set \code{probes = 0} to include all features. Note that providing a numeric vector
#'  for the \code{probes} argument will have \code{plGrid} search across multiple
#'  top features. However, by providing a list of numeric vectors as the \code{probes}
#'  argument, the user can force the default handling of numeric vectors.
#' @param how A character string. Specifies the \code{\link{build}} method to iterate.
#' @param fold A numeric scalar. Specifies the number of folds for cross-validation.
#'  Set \code{fold = 0} to perform leave-one-out cross-validation. Argument passed
#'  to \code{\link{plCV}}. Set \code{fold = NULL} to skip cross-validation altogether.
#' @param aucSkip A logical scalar. Argument passed to \code{\link{calcStats}}.
#' @param verbose A logical scalar. Argument passed to \code{\link{exprso-predict}}.
#' @param ... Arguments passed to the \code{how} method. Unlike the \code{how} wrapper,
#'  \code{plGrid} will accept multiple terms for each argument, supplied as a vector.
#'  \code{plGrid} will train and deploy a classifier for each combination of
#'  parameters provided in this way.
#' @return An \code{\link{ExprsPipeline-class}} object.
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
#' arrays <- splitSample(array, percent.include = 67)
#' array.train <- fsStats(arrays[[1]], probes = 0, how = "t.test")
#' pl <- plGrid(array.train, array.valid = arrays[[2]], how = "buildSVM",
#'              kernel = c("linear", "radial"), cost = 10^(-3:3), gamma = 10^(-3:3))
#' }
#' @export
plGrid <- function(array.train, array.valid = NULL, probes, how, fold = 10,
                   aucSkip = FALSE, verbose = TRUE, ...){

  args <- as.list(substitute(list(...)))[-1]

  if(is.numeric(probes)){

    if(any(probes > nrow(array.train@exprs))){

      message("At least one 'probes' index is too large. Using all probes instead.")
      probes[probes > nrow(array.train@exprs)] <- nrow(array.train@exprs)
    }

    probes <- unique(probes)

  }else if(is.character(probes)){

    # Turn character vector into single list entry
    probes <- as.list(probes)

  } # if list, do nothing here

  # Build grid
  grid <- expand.grid(append(list("probes" = probes), lapply(args, eval)), stringsAsFactors = FALSE)

  # Refine grid for buildSVM
  if(how == "buildSVM"){

    if(!"kernel" %in% names(args)) stop("Uh oh! 'kernel' argument missing!")
    if(!"cost" %in% names(args)) stop("Uh oh! 'cost' argument missing!")
    if("radial" %in% eval(args$kernel)){

      if(!"gamma" %in% names(args)) stop("Uh oh! 'gamma' argument missing!")
      grid[grid$kernel %in% "linear", "gamma"] <- NA
    }

    if("polynomial" %in% eval(args$kernel)){

      if(!"degree" %in% names(args)) stop("Uh oh! 'degree' argument missing!")
      if(!"coef0" %in% names(args)) stop("Uh oh! 'coef0' argument missing!")
      grid[grid$kernel %in% c("linear", "kernel"), "degree"] <- NA
      grid[grid$kernel %in% c("linear", "kernel"), "coef0"] <- NA
    }

    grid <- unique(grid)
  }

  # For each gridpoint in grid
  statistics <- vector("list", nrow(grid))
  models <- vector("list", nrow(grid))
  for(i in 1:nrow(grid)){

    cat("Now building machine at gridpoint:\n")
    print(grid[i, , drop = FALSE])

    # Format gridpoint args to pass along to build do.call
    args <- append(list("object" = array.train), as.list(grid[i, , drop = FALSE]))

    # Build and save model
    args <- lapply(args, unlist)
    model <- do.call(what = how, args = args[!is.na(args)])
    models[[i]] <- model

    # Predict class labels using the provided training set and calculate accuracy
    pred.train <- predict(model, array.train, verbose = verbose)
    stats <- calcStats(pred.train, aucSkip = aucSkip, plotSkip = TRUE)
    colnames(stats) <- paste0("train.", colnames(stats))
    acc <- stats

    # If a validation set is provided
    if(!is.null(array.valid)){

      # Predict class labels using the provided validation set and calculate accuracy
      pred.valid <- predict(model, array.valid, verbose = verbose)
      stats <- calcStats(pred.valid, aucSkip = aucSkip, plotSkip = TRUE)
      colnames(stats) <- paste0("valid.", colnames(stats))
      acc <- data.frame(acc, stats)
    }

    # If 'fold' argument is provided
    if(!is.null(fold)){

      # Perform leave-one-out or v-fold cross-validation
      args <- append(list("how" = how, "fold" = fold), args)
      names(args)[names(args) == "object"] <- "array"
      cv <- do.call(what = plCV, args = args[!is.na(args)])
      acc <- data.frame("fold" = fold, "train.plCV" = cv, acc)
    }

    # Save summary statistics
    statistics[[i]] <- data.frame(grid[i, , drop = FALSE], acc)
  }

  pl <- new("ExprsPipeline",
            summary = do.call(rbind, statistics),
            machs = models
  )

  return(pl)
}
