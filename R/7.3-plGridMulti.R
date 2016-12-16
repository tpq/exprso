###########################################################
### High-throughput classification

#' Perform High-Throughput Classification
#'
#' Trains and deploys multiple classifiers across a vast parameter search space.
#'
#' Unlike \code{plGrid}, the \code{plGridMulti} function accepts a \code{ctrlFS}
#'  argument, allowing for 1-vs-all classification with implicit feature selection.
#   In other words, instead of performing multi-class feature selection followed by
#'  1-vs-all classification, this function divides the data into 1-vs-all bins,
#'  performs a 1-vs-all feature selection for each bin, and then performs a 1-vs-all
#'  classification for that same bin. As such, each \code{ExprsMachine} within the
#'  \code{ExprsModule} will have its own unique feature selection history.
#'
#' Take note, that \code{plGridMulti} does not have built-in \code{plCV} support.
#'  To use \code{plGridMulti} with cross-validation, use \code{plNested}.
#'
#' @inheritParams plGrid
#' @param array.train Specifies the \code{ExprsMulti} object to use as training set.
#' @param array.valid Specifies the \code{ExprsMulti} object to use as validation set.
#' @param ctrlFS A list of arguments handled by \code{\link{ctrlFeatureSelect}}.
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
#' @export
plGridMulti <- function(array.train, array.valid = NULL, ctrlFS, top, how,
                        aucSkip = FALSE, verbose = TRUE, ...){

  if(class(array.train) != "ExprsMulti"){

    stop("This grid search is for multi-class datasets only. Use 'plGrid' instead.")
  }

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
    args <- as.list(grid[i, , drop = FALSE])
    args <- args[!is.na(args)]

    # Build and save N binary models
    multi <- vector("list", length(levels(array.train@annot$defineCase)))
    for(j in 1:length(levels(array.train@annot$defineCase))){

      # If the j-th ExprsMachine would not have any representative cases
      if(all(as.numeric(array.train@annot$defineCase) != j)){

        cat("Missing class ", j, ". Using a NA placeholder instead.\n", sep = "")
        multi[[j]] <- NA

      }else{

        # Turn the ExprsMulti object into the j-th ExprsBinary object
        temp <- array.train
        temp@annot$defineCase <- ifelse(as.numeric(temp$defineCase) == j, "Case", "Control")
        class(temp) <- "ExprsBinary"

        cat("Performing one-vs-all feature selection with class", j, "set as \"Case\".\n")
        if(!"list" %in% lapply(ctrlFS, class)) ctrlFS <- list(ctrlFS)
        for(k in 1:length(ctrlFS)){

          func <- ctrlFS[[k]]$func
          args.fs <- append(list("object" = temp), ctrlFS[[k]][!ctrlFS[[k]] %in% func])
          temp <- do.call(what = func, args = args.fs)
        }

        cat("Performing one-vs-all classification with class", j, "set as \"Case\".\n")
        args.build <- append(list("object" = temp), args)
        multi[[j]] <- do.call(how, args.build)
      }

      model <- new("ExprsModule",
                   preFilter = append(array.train@preFilter, list(top)),
                   reductionModel = append(array.train@reductionModel, list(NA)),
                   mach = multi)
    }

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

    # Save summary statistics
    statistics[[i]] <- data.frame("build" = how, grid[i, , drop = FALSE], acc)
  }

  pl <- new("ExprsPipeline",
            summary = do.call(rbind, statistics),
            machs = models
  )

  return(pl)
}
