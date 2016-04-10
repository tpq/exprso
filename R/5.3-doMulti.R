###########################################################
### doMulti for "1 vs. all" tasks

#' Perform "1 vs. all" Task
#'
#' A function to execute multiple "1 vs. all" binary tasks.
#'
#' \code{doMulti} depends on the total number of levels in the
#'  \code{$defineCase} factor. If a training set is missing any
#'  one of the factor levels (e.g. owing to random cuts during
#'  cross-validation), the \code{ExprsModule} component that
#'  would refer to that class label gets replaced with a NULL
#'  placeholder. This NULL placeholder gets handled as a
#'  special case when predicting with an \code{ExprsModule}.
#'
#' During \code{ExprsModule} class prediction, the absence
#'  of a class during training (i.e. a NULL placeholder)
#'  will prevent an \code{ExprsModule} object from possibly
#'  predicting that class. Instead, it will make a prediction
#'  for class labels included in the training set. However,
#'  classes present in the validation set but missing in the
#'  training set will still undergo prediction and therefore
#'  impact metrics of classifier performance.
#'
#' An \code{ExprsModule} object can only make predictions on
#'  an \code{ExprsMulti} object with the same number of recorded
#'  class labels (i.e. the total number of levels in the
#'  \code{$defineCase} factor). As with all functions included
#'  in this package, all ties get broken using probability
#'  weights proportional to the relative class frequencies
#'  in the training set.
#'
#' @seealso
#' * \code{\link{fs}}, \code{\link{build}}, \code{\link{doMulti}},
#'  \code{\link{reRank}}, \code{\link{exprso-predict}}
#' * \code{\link{plCV}}, \code{\link{plGrid}},
#'  \code{\link{plMonteCarlo}}, \code{\link{plNested}}
#'
#' @export
setGeneric("doMulti",
           function(object, ...) standardGeneric("doMulti")
)

#' @describeIn doMulti Method to execute multiple "1 vs. all" binary tasks.
#' @inheritParams fs
#' @param what A character string. The \code{ExprsBinary} method to execute multiple times.
#' @export
setMethod("doMulti", "ExprsMulti",
          function(object, probes, what, ...){

            args <- as.list(substitute(list(...)))[-1]

            # Initialize multi container
            multi <- vector("list", length(levels(object@annot$defineCase)))

            # Perform N binary tasks
            for(i in 1:length(levels(object@annot$defineCase))){

              # If the i-th ExprsMachine would not have any representative cases
              if(all(!as.numeric(object@annot$defineCase) == i)){

                cat(
                  paste0("Missing a representative of class ", i,
                         ". Object replaced with NULL placeholder.\n")
                )

                multi[[i]] <- NULL

              }else{

                # Turn the ExprsMulti object into the i-th ExprsBinary object
                temp <- object
                temp$defineCase <- ifelse(as.numeric(temp$defineCase) == i, "Case", "Control")
                class(temp) <- "ExprsBinary"

                # Perform the binary task
                cat("Performing a one-vs-all binary task with class", i, "set as \"Case\".\n")
                args.i <- append(list("object" = temp, "probes" = probes), args)
                multi[[i]] <- do.call(what = what, args = args.i)
              }
            }

            return(multi)
          }
)

###########################################################
### reRank to serialize "1 vs. all" feature selection

#' Serialize "1 vs. all" Feature Selection
#'
#' A function to convert multiple feature rank lists from "1 vs. all"
#'  binary feature selection into a single feature rank list.
#'
#' After passing a feature selection method through \code{doMulti},
#'  a set of ranked features gets returned for each one of the
#'  total number of levels in the \code{$defineCase} factor. In
#'  order to proceed with model deployment (at least in the setting
#'  of a conventional pipeline where feature selection occurs
#'  prior to classifier construction), these multiple feature
#'  rankings must get serialized into a single feature rank list.
#'  \code{reRank} accomplishes this by calculating the rank sum
#'  for each feature across all "1 vs. all" feature selection
#'  tasks. Features found in one rank list but not in another
#'  receive a numeric rank equal to one more than the maximum rank
#'  in that feature rank list. The presence of a NULL placeholder
#'  (see: \code{\link{doMulti}}) will not impact \code{reRank}.
#'
#' We note here, however, that a better approach would perform
#'  "1 vs. all" feature selection and classifier construction
#'  simultaneously, rather than "1 vs. all" feature selection
#'  followed by "1 vs. all" classifier construction. We hope to
#'  provide a \code{pl} method for this in the future.
#'
#' We also note here that applying \code{reRank} to data that have
#'  undergone dimension reduction does not make sense.
#'
#' @seealso
#' * \code{\link{fs}}, \code{\link{build}}, \code{\link{doMulti}},
#'  \code{\link{reRank}}, \code{\link{exprso-predict}}
#' * \code{\link{plCV}}, \code{\link{plGrid}},
#'  \code{\link{plMonteCarlo}}, \code{\link{plNested}}
#'
#' @export
reRank <- function(fss){

  # Remove any NULL placeholders
  fss <- fss[!sapply(fss, is.null)]

  # Line up rankings for each doMulti result
  i <- 1
  while(i < length(fss)){

    if(i == 1){

      rank <- data.frame("probe" = rownames(fss[[i]]@exprs),
                         "rank" = 1:nrow(fss[[i]]@exprs))
    }

    rank.next <- data.frame("probe" = rownames(fss[[i + 1]]@exprs),
                            "rank" = 1:nrow(fss[[i + 1]]@exprs))

    rank <- merge(rank, rank.next,
                  by.x = "probe",
                  by.y = "probe",
                  all = TRUE)

    i <- i + 1
  }

  # Clean up rank table
  rownames(rank) <- rank[, "probe"]
  rank <- rank[, !colnames(rank) %in% "probe"]

  # For each class, replace any NAs with 1 more than the maximum rank
  for(col in 1:ncol(rank)){

    maximum <- max(rank[!is.na(rank[, col]), col])
    rank[is.na(rank[, col]), col] <- maximum + 1
  }

  # Add per-class ranks to make final rank list
  final <- rowSums(rank)
  final <- names(final)[order(final)]
}

###########################################################
### Perform ExprsMulti feature selection

#' @describeIn fs Method to perform random feature selection for
#'  multi-class classification.
#' @export
setMethod("fsSample", "ExprsMulti",
          function(object, probes, ...){

            # Call doMulti and make single rank list
            fss <- doMulti(object, probes, what = "fsSample", ...)
            final <- reRank(fss)

            new("ExprsMulti",
                exprs = object@exprs[final,],
                annot = object@annot,
                preFilter = append(object@preFilter, list(final)),
                reductionModel = append(object@reductionModel, list(NA)))
          }
)

#' @describeIn fs Method to perform statistics based feature
#'  selection for multi-class classification.
#' @export
setMethod("fsStats", "ExprsMulti",
          function(object, probes, ...){

            # Call doMulti and make single rank list
            fss <- doMulti(object, probes, what = "fsStats", ...)
            final <- reRank(fss)

            new("ExprsMulti",
                exprs = object@exprs[final,],
                annot = object@annot,
                preFilter = append(object@preFilter, list(final)),
                reductionModel = append(object@reductionModel, list(NA)))
          }
)

###########################################################
### Perform ExprsMulti build

#' @describeIn build Method to build classifiers for
#'  multi-class classification.
#' @export
setMethod("buildNB", "ExprsMulti",
          function(object, probes, ...){

            # Pass arguments to doMulti
            machs <- doMulti(object, probes, what = "buildNB", ...)

            # Carry through and append fs history as stored in the ExprsArray object
            new("ExprsModule",
                preFilter = append(object@preFilter, list(probes)),
                reductionModel = append(object@reductionModel, list(NA)),
                mach = machs)
          }
)

#' @describeIn build Method to build classifiers for
#'  multi-class classification.
#' @export
setMethod("buildLDA", "ExprsMulti",
          function(object, probes, ...){

            # Pass arguments to doMulti
            machs <- doMulti(object, probes, what = "buildLDA", ...)

            # Carry through and append fs history as stored in the ExprsArray object
            new("ExprsModule",
                preFilter = append(object@preFilter, list(probes)),
                reductionModel = append(object@reductionModel, list(NA)),
                mach = machs)
          }
)

#' @describeIn build Method to build classifiers for
#'  multi-class classification.
#' @export
setMethod("buildSVM", "ExprsMulti",
          function(object, probes, ...){

            # Pass arguments to doMulti
            machs <- doMulti(object, probes, what = "buildSVM", ...)

            # Carry through and append fs history as stored in the ExprsArray object
            new("ExprsModule",
                preFilter = append(object@preFilter, list(probes)),
                reductionModel = append(object@reductionModel, list(NA)),
                mach = machs)
          }
)

#' @describeIn build Method to build classifiers for
#'  multi-class classification.
#' @export
setMethod("buildANN", "ExprsMulti",
          function(object, probes, ...){

            # Pass arguments to doMulti
            machs <- doMulti(object, probes, what = "buildANN", ...)

            # Carry through and append fs history as stored in the ExprsArray object
            new("ExprsModule",
                preFilter = append(object@preFilter, list(probes)),
                reductionModel = append(object@reductionModel, list(NA)),
                mach = machs)
          }
)

#' @describeIn build Method to build classifiers for
#'  multi-class classification.
#' @export
setMethod("buildRF", "ExprsMulti",
          function(object, probes, ...){

            # Pass arguments to doMulti
            machs <- doMulti(object, probes, what = "buildRF", ...)

            # Carry through and append fs history as stored in the ExprsArray object
            new("ExprsModule",
                preFilter = append(object@preFilter, list(probes)),
                reductionModel = append(object@reductionModel, list(NA)),
                mach = machs)
          }
)
