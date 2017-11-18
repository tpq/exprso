###########################################################
### High-throughput classification

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

    # Turn character vector into single list entry
    top <- as.list(top)

  } # if list, do nothing here

  args <- getArgs(...)
  grid <- expand.grid(append(list("top" = top), lapply(args, eval)), stringsAsFactors = FALSE)

  # Refine grid for buildSVM
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

#' Perform High-Throughput Classification
#'
#' Trains and deploys multiple classifiers across a vast parameter search space.
#'
#' \code{plGrid} will \code{\link{build}} and \code{\link{exprso-predict}} for
#'  each combination of parameters provided as additional arguments (\code{...}).
#'  When using \code{plGrid}, supplying a numeric vector as the \code{top}
#'  argument will train and deploy a classifier of each mentioned size for
#'  each combination of parameters provided in \code{...}. To skip validation set
#'  prediction, use \code{array.valid = NULL}. Either way, this function returns an
#'  \code{\link{ExprsPipeline-class}} object which contains the build parameters,
#'  classifier performances, and \code{ExprsModel} objects for all trained models.
#'
#' \code{plGrid} will perform v-fold or leave-one-out cross-validation for the
#'  using \code{\link{plCV}}. The argument \code{fold} specifies the number
#'  of v-folds to use during cross-validation. Set \code{fold = 0} to perform
#'  leave-one-out cross-validation. This approach to cross-validation
#'  will work for \code{ExprsBinary} and \code{ExprsMulti} objects alike. The
#'  performance metric used to measure cross-validation accuracy is the
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
#' @param top A numeric scalar or character vector. A numeric scalar indicates
#'  the number of top features that should undergo feature selection. A character vector
#'  indicates specifically which features by name should undergo feature selection.
#'  Set \code{top = 0} to include all features. Note that providing a numeric vector
#'  for the \code{top} argument will have \code{plGrid} search across multiple
#'  top features. However, by providing a list of numeric vectors as the \code{top}
#'  argument, the user can force the default handling of numeric vectors.
#' @param how A character string. Specifies the \code{\link{build}} method to iterate.
#' @param fold A numeric scalar. Specifies the number of folds for cross-validation.
#'  Set \code{fold = 0} to perform leave-one-out cross-validation. Argument passed
#'  to \code{\link{plCV}}. Set \code{fold = NULL} to skip cross-validation altogether.
#' @param aucSkip A logical scalar. Argument passed to \code{\link{calcStats}}.
#' @param verbose A logical scalar. Argument passed to \code{\link{exprso-predict}}.
#' @param ... Arguments passed to the \code{how} method. Unlike the \code{build} method,
#'  \code{plGrid} allows multiple parameters for each argument, supplied as a vector.
#'  See Details.
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
#' arrays <- splitSample(array, percent.include = 67)
#' array.train <- fsStats(arrays[[1]], top = 0, how = "t.test")
#' pl <- plGrid(array.train, array.valid = arrays[[2]], how = "buildSVM",
#'              kernel = c("linear", "radial"), cost = 10^(-3:3), gamma = 10^(-3:3))
#' }
#' @export
plGrid <- function(array.train, array.valid = NULL, top, how, fold = 10,
                   aucSkip = FALSE, verbose = TRUE, ...){

  if(missing(how)){

    stop("Uh oh! You must provide a valid build method for the 'how' argument.")
  }

  # For each gridpoint in grid
  grid <- makeGridFromArgs(array.train = array.train, top = top, how = how, ...)
  grid <- grid[, !colnames(grid) %in% "plotSkip", drop = FALSE]
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
    statistics[[i]] <- data.frame("build" = how, grid[i, , drop = FALSE], acc)
  }

  pl <- new("ExprsPipeline",
            summary = do.call(rbind, statistics),
            machs = models
  )

  return(pl)
}
