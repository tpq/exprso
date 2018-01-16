#' Perform High-Throughput Machine Learning
#'
#' Trains and deploys models across a vast parameter search space.
#'
#' \code{plGrid} will \code{\link{build}} and \code{\link{exprso-predict}} for
#'  each combination of parameters provided as additional arguments (\code{...}).
#'  When using \code{plGrid}, supplying a numeric vector as the \code{top}
#'  argument will train and deploy a model of each mentioned size for
#'  each combination of parameters provided in \code{...}. To skip validation set
#'  prediction, use \code{array.valid = NULL}. Either way, this function returns an
#'  \code{\link{ExprsPipeline-class}} object which contains a summary of the build
#'  parameters and the models themselves. The argument \code{fold} controls
#'  cross-validation via \code{\link{plCV}}.
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
#' @return An \code{\link{ExprsPipeline-class}} object.
#' @export
plGrid <- function(array.train, array.valid = NULL, top, how, fold = 10,
                   aucSkip = FALSE, verbose = FALSE, ...){

  if(missing(how)){

    stop("Uh oh! You must provide a valid build method for the 'how' argument.")
  }

  # For each gridpoint in grid
  grid <- makeGridFromArgs(array.train = array.train, top = top, how = how, ...)
  grid <- grid[, !colnames(grid) %in% "plotSkip", drop = FALSE]
  statistics <- vector("list", nrow(grid))
  models <- vector("list", nrow(grid))
  for(i in 1:nrow(grid)){

    if(verbose){
      cat("Now building machine at gridpoint:\n")
      print(grid[i, , drop = FALSE])
    }

    # Format gridpoint args to pass along to build do.call
    args <- append(list("object" = array.train), as.list(grid[i, , drop = FALSE]))

    # Build and save model
    args <- lapply(args, unlist)
    model <- do.call(what = how, args = args[!is.na(args)])
    models[[i]] <- model

    # Predict class labels using the provided training set and calculate accuracy
    pred.train <- predict(model, array.train, verbose = verbose)
    stats <- calcStats(pred.train, aucSkip = aucSkip, plotSkip = TRUE, verbose = FALSE)
    colnames(stats) <- paste0("train.", colnames(stats))
    acc <- stats

    # If a validation set is provided
    if(!is.null(array.valid)){

      # Predict class labels using the provided validation set and calculate accuracy
      pred.valid <- predict(model, array.valid, verbose = verbose)
      stats <- calcStats(pred.valid, aucSkip = aucSkip, plotSkip = TRUE, verbose = FALSE)
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
