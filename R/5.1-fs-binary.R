###########################################################
### Define generic functions

#' @name fs
#' @rdname fs
#'
#' @title Perform Feature Selection
#'
#' @description A collection of functions to select features.
#'
#' @details
#'
#' Considering the high-dimensionality of most genomic datasets, it is prudent and often necessary
#'  to prioritize which features to include during classifier construction. There exists a myriad of
#'  ways to perform the task of feature selection. This package provides methods for some of the most
#'  frequently used feature selection methods. Each function works as a self-contained wrapper that
#'  (1) pre-processes the ExprsArray input, (2) performs the feature selection, and (3) returns an
#'  ExprsArray output with an updated feature selection history. The user may deploy, in tandem, any
#'  number of these functions in whatever order they choose, limited only by computational power and
#'  imagination. These feature selection histories get passed along at every step of the way until
#'  they eventually get used in order to pre-process an unlabelled dataset during classifier deployment
#'  (i.e. prediction). In the spirit of open source programming, we encourage users to submit their
#'  own feature selection functions, modeled after those provided in this library.
#'
#' All feature selection methods here perform feature selection prior to classifier construction. For
#'  \code{ExprsBinary} methods, this is typically the desired approach. However, the exprso package
#'  also provides \code{pl} methods which automate integrated machine learning pipelines.
#'  However, these feature selection methods do not necessarily generalize to multi-class classification.
#'  As such, the \code{ExprsMulti} methods instead harness the \code{\link{doMulti}} function
#'  to perform "1 vs. all" binary feature selection, aggregating the final results with
#'  \code{\link{reRank}}. A better approach would perform "1 vs. all" feature selection and classifier
#'  construction simultaneously, rather than "1 vs. all" feature selection followed by "1 vs. all"
#'  classifier construction. We hope to provide a \code{pl} method for this in the future.
#'
#' For all feature selection methods, \code{@@preFilter} and \code{@@reductionModel} stores the
#'  feature selection and dimension reduction history, respectively. This history gets passed
#'  along to prepare the test or validation set during model deployment, ensuring that these
#'  sets undergo the same feature selection and dimension reduction in the appropriate steps.
#'  Under the scenarios where users plan to apply multiple feature selection or dimension
#'  reduction steps, the \code{probes} argument manages which features (e.g. gene probes) to
#'  send through each feature selection or dimension reduction procedure. For \code{probes},
#'  a numeric scalar indicates the number of top features to use, while a character vector
#'  indicates specifically which features to use. In this way, the user sets which features
#'  to feed INTO the \code{fs} method (NOT which features the user expects OUT). The example
#'  below shows how to apply dimension reduction to the top 50 features as selected by the
#'  Student's t-test. Set \code{probes = 0} to pass all features through an \code{fs} method.
#'
#' Note that \code{fsMrmre} crashes when supplied a very large \code{feature_count} argument
#'  owing to its implementation in the imported package \code{mRMRe}.
#'
#' @param object Specifies the \code{ExprsArray} object to undergo feature selection.
#' @param probes A numeric scalar or character vector. A numeric scalar indicates
#'  the number of top features that should undergo feature selection. A character vector
#'  indicates specifically which features by name should undergo feature selection.
#'  Set \code{probes = 0} to include all features. A numeric vector can also be used
#'  to indicate specific features by location, similar to a character vector.
#' @param how Specifics which function to call in \code{fsStats}. Recognized arguments
#'  include \code{"t.test"} and \code{"ks.test"}.
#' @param ... Arguments passed to the respective wrapped function.
#' @return Returns an \code{ExprsArray} object.
#'
#' @seealso
#' \code{\link{fs}}\cr
#' \code{\link{build}}\cr
#' \code{\link{doMulti}}\cr
#' \code{\link{reRank}}\cr
#' \code{\link{exprso-predict}}\cr
#' \code{\link{plCV}}\cr
#' \code{\link{plGrid}}\cr
#' \code{\link{plMonteCarlo}}\cr
#' \code{\link{plNested}}
#'
#' @examples
#' \dontrun{
#' library(golubEsets)
#' data(Golub_Merge)
#' array <- arrayEset(Golub_Merge, colBy = "ALL.AML", include = list("ALL", "AML"))
#' array <- modFilter(array, 20, 16000, 500, 5) # pre-filter Golub ala Deb 2003
#' array <- modTransform(array) # lg transform
#' array <- modNormalize(array, c(1, 2)) # normalize gene and subject vectors
#' arrays <- splitSample(array, percent.include = 67)
#' array.train <- fsStats(arrays[[1]], probes = 0, how = "t.test")
#' array.train <- fsPrcomp(array.train, probes = 50)
#' mach <- buildSVM(array.train, probes = 5, kernel = "linear", cost = 1)
#' }
NULL

#' @rdname fs
#' @export
setGeneric("fsSample",
           function(object, probes = 0, ...) standardGeneric("fsSample")
)

#' @rdname fs
#' @export
setGeneric("fsStats",
           function(object, probes = 0, ...) standardGeneric("fsStats")
)

#' @rdname fs
#' @export
setGeneric("fsPrcomp",
           function(object, probes = 0, ...) standardGeneric("fsPrcomp")
)

#' @rdname fs
#' @export
setGeneric("fsPathClassRFE",
           function(object, probes = 0, ...) standardGeneric("fsPathClassRFE")
)

#' @rdname fs
#' @export
setGeneric("fsEbayes",
           function(object, probes = 0, ...) standardGeneric("fsEbayes")
)

#' @rdname fs
#' @export
setGeneric("fsMrmre",
           function(object, probes = 0, ...) standardGeneric("fsMrmre")
)

###########################################################
### Select features

#' Workhorse for fs Methods
#'
#' Used as a back-end wrapper for creating new fs methods.
#'
#' If the uniqueFx returns a character vector, it is assumed
#'  that the fs method is a feature selection only. If the
#'  uniqueFx returns a list, it is assumed that the fs method
#'  is a reduction model method only.
#'
#' @inheritParams fs
#' @param uniqueFx A function call unique to that fs method.
#' @return Returns an \code{ExprsArray} object.
#'
#' @export
fs. <- function(object, probes, uniqueFx, ...){

  if(class(probes) == "numeric"){

    if(length(probes) == 1){

      if(probes > nrow(object@exprs)) probes <- 0
      if(probes == 0) probes <- nrow(object@exprs)
      probes <- rownames(object@exprs[1:probes, ])

    }else{

      probes <- rownames(object@exprs[probes, ])
    }
  }

  data <- t(object@exprs[probes, ])
  final <- do.call("uniqueFx", list(data, probes, ...))

  if(class(final) == "character"){

    array <- new(class(object), exprs = object@exprs[final,], annot = object@annot,
                 preFilter = append(object@preFilter, list(final)),
                 reductionModel = append(object@reductionModel, list(NA))
    )

  }else if(class(final) == "list"){

    array <- new(class(object), exprs = final[[1]], annot = object@annot,
                 preFilter = append(object@preFilter, list(probes)),
                 reductionModel = append(object@reductionModel, list(final[[2]]))
    )

  }else{

    stop("Uh oh! DEBUG ERROR: 002")
  }
}

#' @rdname fs
#' @section Methods (by generic):
#' \code{fsSample:} Method to perform random feature selection using base::sample.
#' @export
setMethod("fsSample", "ExprsBinary",
          function(object, probes, ...){ #args to sample

            fs.(object, probes,
                uniqueFx = function(data, probes, ...){

                  sample(probes, ...)
                }, ...)
          }
)

#' @rdname fs
#' @section Methods (by generic):
#' \code{fsStats:} Method to perform statistics based feature selection using stats::t.test and others.
#' @importFrom stats t.test ks.test
#' @export
setMethod("fsStats", "ExprsBinary",
          function(object, probes, how = "t.test", ...){ # args to t.test or ks.test

            if(how == "t.test"){

              fs.(object, probes,
                  uniqueFx = function(data, probes, ...){

                    # Prepare data for statistical tests
                    cases <- object@annot$defineCase %in% "Case"
                    conts <- object@annot$defineCase %in% "Control"
                    p <- vector("numeric", length(probes))

                    for(i in 1:length(probes)){

                      tryCatch(
                        {
                          p[i] <- t.test(object@exprs[probes[i], cases],
                                         object@exprs[probes[i], conts], ...)$p.value

                        }, error = function(e){

                          cat("fsStats failed for feature: ", probes[i], ". Setting p(x)=1...\n")
                          p[i] <- 1
                        })
                    }

                    probes[order(p)]
                  }, ...)

            }else if(how == "ks.test"){

              fs.(object, probes,
                  uniqueFx = function(data, probes, ...){

                    # Prepare data for statistical tests
                    cases <- object@annot$defineCase %in% "Case"
                    conts <- object@annot$defineCase %in% "Control"
                    p <- vector("numeric", length(probes))

                    for(i in 1:length(probes)){

                      tryCatch(
                        {
                          p[i] <- ks.test(object@exprs[probes[i], cases],
                                          object@exprs[probes[i], conts], ...)$p.value

                        }, error = function(e){

                          cat("fsStats failed for feature: ", probes[i], ". Setting p(x)=1...\n")
                          p[i] <- 1
                        })
                    }

                    probes[order(p)]
                  }, ...)

            }else{

              stop("Uh oh! Provided 'how' argument not recognized!")
            }
          }
)

#' @rdname fs
#' @section Methods (by generic):
#' \code{fsPrcomp:} Method to perform principal components analysis using stats::prcomp.
#' @importFrom stats prcomp
#' @export
setMethod("fsPrcomp", "ExprsBinary",
          function(object, probes, ...){ # args to prcomp

            fs.(object, probes,
                uniqueFx = function(data, probes, ...){

                  # ATTENTION: We want dependent variables as columns
                  # NOTE: as.data.frame will not rename columns
                  #  -- I don't understand this note? 19/09/16
                  reductionModel <- prcomp(data.frame(data), ...)

                  # ATTENTION: The value of predict(reductionModel, data) equals $x
                  # This information will automatically distill the data
                  #  when calling modHistory
                  list(t(reductionModel$x),
                       reductionModel)
                }, ...)
          }
)

#' @rdname fs
#' @section Methods (by generic):
#' \code{fsPathClassRFE:} Method to perform SVM-RFE feature selection using pathClass::fit.rfe.
#' @importFrom pathClass fit.rfe
#' @export
setMethod("fsPathClassRFE", "ExprsBinary",
          function(object, probes, ...){ # args to fit.rfe

            fs.(object, probes,
                uniqueFx = function(data, probes, ...){

                  if(!requireNamespace("kernlab", quietly = TRUE)){
                    stop("Uh oh! This fs method depends on kernlab! ",
                         "Try running: install.packages('kernlab')")

                  }else{

                    # Removing this will break the method:
                    library(kernlab)
                  }

                  # Set up "make.names" key for improper @exprs row.names
                  labels <- factor(object@annot[rownames(data), "defineCase"], levels = c("Control", "Case"))
                  key <- data.frame("old" = colnames(data), "new" = make.names(colnames(data)))

                  # NOTE: RFE as assembled by pathClass is via a linear kernel only
                  # NOTE: By default, fit.rfe iterates through C = 10^c(-3:3)
                  # Run fit.rfe()
                  rfe <- pathClass::fit.rfe(x = data, y = labels, ...)

                  # Use "make.names" key to return to original row.names
                  final <- merge(data.frame("new" = rfe$features), key, sort = FALSE)$old
                  if(length(final) < 2) stop("Uh oh! fsPathClassRFE did not find enough features!")
                  as.character(final)
                }, ...)
          }
)

#' @rdname fs
#' @section Methods (by generic):
#' \code{fsEbayes:} Method to perform empiric Bayes feature selection using limma::ebayes.
#' @importFrom limma ebayes lmFit
#' @export
setMethod("fsEbayes", "ExprsBinary",
          function(object, probes, ...){ # args to ebayes

            fs.(object, probes,
                uniqueFx = function(data, probes, ...){

                  design <- as.matrix(ifelse(object@annot$defineCase == "Case", 1, 0))
                  colnames(design) <- "CaseVCont"
                  fit <- limma::lmFit(object@exprs[probes, ], design)
                  ebaye <- limma::ebayes(fit, ...)
                  rownames(ebaye$p.value)[order(ebaye$p.value[, 1])]
                }, ...)
          }
)

#' @rdname fs
#' @section Methods (by generic):
#' \code{fsMrme:} Method to perform mRMR feature selection using mRMRe::mRMR.classic.
#' @importFrom mRMRe mRMR.classic mRMR.data solutions featureNames
#' @export
setMethod("fsMrmre", "ExprsBinary",
          function(object, probes, ...){ # args to mRMR.classic

            fs.(object, probes,
                uniqueFx = function(data, probes, ...){

                  args <- getArgs(...)
                  args <- defaultArg("target_indices", 1, args)
                  args <- defaultArg("feature_count", 64, args)

                  # Set up "make.names" key for improper @exprs row.names
                  key <- data.frame("old" = colnames(data), "new" = make.names(colnames(data)))

                  labels <- as.numeric(object@annot$defineCase == "Case")
                  mRMRdata <- mRMRe::mRMR.data(data = data.frame(labels, data))
                  args <- append(list("data" = mRMRdata), args)
                  mRMRout <- do.call(mRMRe::mRMR.classic, args)

                  # Sort probes
                  final <- as.vector(
                    apply(mRMRe::solutions(mRMRout)[[1]], 2, function(x, y){ return(y[x])},
                          y = mRMRe::featureNames(mRMRdata))
                  )

                  # Use "make.names" key to return to original row.names
                  final <- merge(data.frame("new" = final), key, sort = FALSE)$old
                  as.character(final)
                }, ...)
          }
)
