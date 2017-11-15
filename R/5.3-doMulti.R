###########################################################
### doMulti for "1 vs. all" tasks

#' Perform "1 vs. all" Task
#'
#' A function to execute multiple "1 vs. all" binary tasks.
#'
#' \code{doMulti} depends on the total number of levels in the
#'  \code{$defineCase} factor. If a training set is missing any
#'  one of the factor levels (e.g., owing to random cuts during
#'  cross-validation), the \code{ExprsModule} component that
#'  would refer to that class label gets replaced with an NA
#'  placeholder. This NA placeholder gets handled as a
#'  special case when predicting with an \code{ExprsModule}.
#'
#' During \code{ExprsModule} class prediction, the absence
#'  of a class during training (i.e., an NA placeholder)
#'  will prevent an \code{ExprsModule} object from possibly
#'  predicting that class in a validation set. Rather, an
#'  \code{ExprsModule} can only make predictions about class
#'  labels that it "knows". However, all "unknown" classes
#'  in the validation set (i.e., those missing from the training
#'  set) still impact metrics of classifier performance.
#'
#' An \code{ExprsModule} object can only make predictions on
#'  an \code{ExprsMulti} object with the same number of recorded
#'  class labels (i.e., the total number of levels in the
#'  \code{$defineCase} factor). As with all functions included
#'  in this package, all ties get broken using probability
#'  weights proportional to the relative class frequencies
#'  in the training set.
#'
#' @inheritParams fs.
#' @param method A character string. The \code{ExprsBinary} method to execute multiple times.
#' @return A list of the results given by \code{method}.
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
setGeneric("doMulti",
           function(object, top, method, ...) standardGeneric("doMulti")
)

#' @describeIn doMulti Method to execute multiple "1 vs. all" binary tasks.
#' @export
setMethod("doMulti", "ExprsMulti",
          function(object, top, method, ...){

            # Perform N binary tasks
            args <- getArgs(...)
            multi <- vector("list", length(levels(object@annot$defineCase)))
            for(i in 1:length(levels(object@annot$defineCase))){

              # If the i-th ExprsMachine would not have any representative cases
              if(all(as.numeric(object@annot$defineCase) != i)){

                cat("Missing class ", i, ". Using a NA placeholder instead.\n", sep = "")
                multi[[i]] <- NA

              }else{

                # Turn the ExprsMulti object into the i-th ExprsBinary object
                temp <- object
                temp@annot$defineCase <- ifelse(as.numeric(temp$defineCase) == i, "Case", "Control")
                class(temp) <- "ExprsBinary"

                # Perform the binary task
                cat("Performing a one-vs-all binary task with class", i, "set as \"Case\".\n")
                args.i <- append(list("object" = temp, "top" = top), args)
                multi[[i]] <- do.call(method, args.i)
              }
            }

            return(multi)
          }
)

###########################################################
### reRank to serialize "1 vs. all" feature selection

#' Serialize "1 vs. all" Feature Selection
#'
#' This experimental function converts multiple feature rank lists,
#'  derived from "1 vs. all" binary feature selection, into a single
#'  feature rank list. This function is not in use in this package.
#'
#' After passing a feature selection method through \code{doMulti},
#'  a set of ranked features gets returned for each one of the
#'  total number of levels in the \code{$defineCase} factor. In
#'  order to proceed with model deployment (at least in the setting
#'  of a conventional pipeline where feature selection occurs
#'  prior to classifier construction), multiple feature rankings
#'  would need to get serialized into a single feature rank list.
#'  \code{reRank} accomplishes this by calculating the rank sum
#'  for each feature across all "1 vs. all" feature selection
#'  tasks. Features found in one rank list but not in another
#'  receive a numeric rank equal to one more than the maximum rank
#'  in that feature rank list. The presence of a NA placeholder
#'  (see: \code{\link{doMulti}}) will not impact \code{reRank}.
#'
#' We note here, however, that a better approach would deploy
#'  "1 vs. all" feature selection and classifier construction
#'  simultaneously, rather than "1 vs. all" feature selection
#'  followed by "1 vs. all" classifier construction. This is
#'  now implemented as \code{\link{plGridMulti}}.
#'
#' @param fss The result of a \code{doMulti} function call.
#' @return A vector of re-ranked features. See Details.
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
reRank <- function(fss){

  # Remove any NULL placeholders
  fss <- fss[!sapply(fss, is.null)]

  # Line up rankings for each doMulti result
  i <- 1
  while(i < length(fss)){

    if(i == 1){

      rank <- data.frame("feat" = rownames(fss[[i]]@exprs),
                         "rank" = 1:nrow(fss[[i]]@exprs))
    }

    rank.next <- data.frame("feat" = rownames(fss[[i + 1]]@exprs),
                            "rank" = 1:nrow(fss[[i + 1]]@exprs))

    rank <- merge(rank, rank.next,
                  by.x = "feat",
                  by.y = "feat",
                  all = TRUE)

    i <- i + 1
  }

  # Clean up rank table
  rownames(rank) <- rank[, "feat"]
  rank <- rank[, !colnames(rank) %in% "feat"]

  # For each class, replace any NAs with 1 more than the maximum rank
  for(col in 1:ncol(rank)){

    maximum <- max(rank[!is.na(rank[, col]), col])
    rank[is.na(rank[, col]), col] <- maximum + 1
  }

  # Add per-class ranks to make final rank list
  final <- rowSums(rank)
  final <- names(final)[order(final)]

  return(final)
}
