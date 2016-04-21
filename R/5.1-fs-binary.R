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
#'  also provides \code{\link{pl}} methods which automate integrated machine learning pipelines.
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
           function(object, ...) standardGeneric("fsSample")
)

#' @rdname fs
#' @export
setGeneric("fsStats",
           function(object, ...) standardGeneric("fsStats")
)

#' @rdname fs
#' @export
setGeneric("fsPrcomp",
           function(object, ...) standardGeneric("fsPrcomp")
)

#' @rdname fs
#' @export
setGeneric("fsPenalizedSVM",
           function(object, ...) standardGeneric("fsPenalizedSVM")
)

#' @rdname fs
#' @export
setGeneric("fsPathClassRFE",
           function(object, ...) standardGeneric("fsPathClassRFE")
)

#' @rdname fs
#' @export
setGeneric("fsEbayes",
           function(object, ...) standardGeneric("fsEbayes")
)

#' @rdname fs
#' @export
setGeneric("fsMrmre",
           function(object, ...) standardGeneric("fsMrmre")
)

###########################################################
### Select features

#' @rdname fs
#' @section Methods (by generic):
#' \code{fsSample:} Method to perform random feature selection using base::sample.
#'
#' @param object Specifies the \code{ExprsArray} object to undergo feature selection.
#' @param probes A numeric scalar or character vector. A numeric scalar indicates
#'  the number of top features that should undergo feature selection. A character vector
#'  indicates specifically which features by name should undergo feature selection.
#'  Set \code{probes = 0} to include all features.
#' @param ... Arguments passed to the respective wrapped function.
#' @return Returns an \code{ExprsArray} object.
#'
#' @export
setMethod("fsSample", "ExprsBinary",
          function(object, probes, ...){ #args to ebayes

            # Convert 'numeric' probe argument to 'character' probe vector
            if(class(probes) == "numeric"){

              if(probes == 0) probes <- nrow(object@exprs)
              probes <- rownames(object@exprs[1:probes, ])
            }

            # Build data using supplied 'character' probe vector
            if(class(probes) == "character"){

              data <- t(object@exprs[probes, ])
            }

            # Randomly sample probes
            final <- sample(probes)

            array <- new("ExprsBinary",
                         exprs = object@exprs[final,],
                         annot = object@annot,
                         preFilter = append(object@preFilter, list(final)),
                         reductionModel = append(object@reductionModel, list(NA))
            )

            return(array)
          }
)

#' @rdname fs
#' @section Methods (by generic):
#' \code{fsStats:} Method to perform statistics based feature selection using stats::t.test and others.
#'
#' @param how Specifics which function to call in \code{fsStats}. Accepted arguments
#'  include \code{"t.test"}, \code{"ks.test"}, and \code{"ks.boot"}.
#'
#' @import Matching
#' @export
setMethod("fsStats", "ExprsBinary",
          function(object, probes, how = "t.test", ...){ # args to ks.test, ks.boot, or t.test

            # Convert 'numeric' probe argument to 'character' probe vector
            if(class(probes) == "numeric"){

              if(probes == 0) probes <- nrow(object@exprs)
              probes <- rownames(object@exprs[1:probes, ])
              data <- t(object@exprs[probes, ])
            }

            # Build data using supplied 'character' probe vector
            if(class(probes) == "character"){

              data <- t(object@exprs[probes, ])
            }

            # Retrieve case and control subjectIDs
            cases <- object@annot$defineCase %in% "Case"
            conts <- object@annot$defineCase %in% "Control"

            # Initialize p-value container
            p <- vector("numeric", length(probes))

            for(i in 1:length(probes)){

              if(how == "ks.test"){

                p[i] <- ks.test(object@exprs[i, cases], object@exprs[i, conts], ...)$p.value

              }else if(how == "ks.boot"){

                p[i] <- Matching::ks.boot(object@exprs[i, cases], object@exprs[i, conts], ...)$ks.boot.pvalue

              }else if(how == "t.test"){

                p[i] <- t.test(object@exprs[i, cases], object@exprs[i, conts], ...)$p.value

              }else{

                stop("Uh oh! Provided 'how' argument not recognized!")}
            }

            final <- probes[order(p)]
            array <- new("ExprsBinary",
                         exprs = object@exprs[final,],
                         annot = object@annot,
                         preFilter = append(object@preFilter, list(final)),
                         reductionModel = append(object@reductionModel, list(NA))
            )

            return(array)
          }
)

#' @rdname fs
#' @section Methods (by generic):
#' \code{fsPrcomp:} Method to perform principal components analysis using stats::prcomp.
#' @export
setMethod("fsPrcomp", "ExprsBinary",
          function(object, probes, ...){ # args to prcomp

            # Convert 'numeric' probe argument to 'character' probe vector
            if(class(probes) == "numeric"){

              if(probes == 0) probes <- nrow(object@exprs)
              probes <- rownames(object@exprs[1:probes, ])
              data <- t(object@exprs[probes, ])
            }

            # Build data using supplied 'character' probe vector
            if(class(probes) == "character"){

              data <- t(object@exprs[probes, ])
            }

            # NOTE: as.data.frame will not rename columns
            data <- data.frame(data)

            # ATTENTION: We want dependent variables as columns
            reductionModel <- prcomp(data, ...)

            cat("\nDimension reduction model summary:\n\n")
            print(summary(reductionModel))

            # ATTENTION: The value of predict(reductionModel, data) equals $x
            # @preFilter stores probes used to build reductionModel (i.e. as passed on by 'probes' argument)
            # This information will automatically distill the data when calling modHistory
            array <- new("ExprsBinary",
                         exprs = t(reductionModel$x),
                         annot = object@annot,
                         preFilter = append(object@preFilter, list(probes)),
                         reductionModel = append(object@reductionModel, list(reductionModel))
            )

            return(array)
          }
)

#' @rdname fs
#' @section Methods (by generic):
#' \code{fsPenalizedSVM:} Method to perform penalizedSVM feature selection using penalizedSVM::svm.fs.
#' @import penalizedSVM
#' @export
setMethod("fsPenalizedSVM", "ExprsBinary",
          function(object, probes, ...){ # args to svm.fs

            # Convert 'numeric' probe argument to 'character' probe vector
            if(class(probes) == "numeric"){

              if(probes == 0) probes <- nrow(object@exprs)
              probes <- rownames(object@exprs[1:probes, ])
              data <- t(object@exprs[probes, ])
            }

            # Build data using supplied 'character' probe vector
            if(class(probes) == "character"){

              data <- t(object@exprs[probes, ])
            }

            # Set labels as factor
            labels <- factor(object@annot[rownames(data), "defineCase"], levels = c("Control", "Case"))
            levels(labels) <- c(-1, 1)

            # Run svm.fs()
            fs <- penalizedSVM::svm.fs(x = data, y = labels, ...)

            final <- names(fs$model$w)

            if(length(final) < 2) stop("Uh oh! fsPenalizedSVM did not find enough features!")
            array <- new("ExprsBinary",
                         exprs = object@exprs[final, ],
                         annot = object@annot,
                         preFilter = append(object@preFilter, list(final)),
                         reductionModel = append(object@reductionModel, list(NA))
            )

            return(array)
          }
)

#' @rdname fs
#' @section Methods (by generic):
#' \code{fsPathClassRFE:} Method to perform SVM-RFE feature selection using pathClass::fit.rfe.
#' @import pathClass
#' @export
setMethod("fsPathClassRFE", "ExprsBinary",
          function(object, probes, ...){ # args to fit.rfe

            # Convert 'numeric' probe argument to 'character' probe vector
            if(class(probes) == "numeric"){

              if(probes == 0) probes <- nrow(object@exprs)
              probes <- rownames(object@exprs[1:probes, ])
              data <- t(object@exprs[probes, ])
            }

            # Build data using supplied 'character' probe vector
            if(class(probes) == "character"){

              data <- t(object@exprs[probes, ])
            }

            # Set labels as factor
            labels <- factor(object@annot[rownames(data), "defineCase"], levels = c("Control", "Case"))

            # Set up "make.names" key for improper @exprs row.names
            key <- data.frame("old" = colnames(data), "new" = make.names(colnames(data)))

            # NOTE: RFE as assembled by pathClass is via a linear kernel only
            # NOTE: By default, fit.rfe iterates through C = 10^c(-3:3)
            # Run fit.rfe()
            rfe <- pathClass::fit.rfe(x = data, y = labels, ...)

            # Sort probes
            final <- rfe$features

            # Use "make.names" key to return to original row.names
            final <- merge(data.frame("new" = final), key)$old

            if(length(final) < 2) stop("Uh oh! fsPathClassRFE did not find enough features!")
            array <- new("ExprsBinary",
                         exprs = object@exprs[final, ],
                         annot = object@annot,
                         preFilter = append(object@preFilter, list(final)),
                         reductionModel = append(object@reductionModel, list(NA))
            )

            return(array)
          }
)

#' @rdname fs
#' @section Methods (by generic):
#' \code{fsEbayes:} Method to perform empiric Bayes feature selection using limma::ebayes.
#' @import limma
#' @export
setMethod("fsEbayes", "ExprsBinary",
          function(object, probes, ...){ # args to ebayes

            # Convert 'numeric' probe argument to 'character' probe vector
            if(class(probes) == "numeric"){

              if(probes == 0) probes <- nrow(object@exprs)
              probes <- rownames(object@exprs[1:probes, ])
            }

            # Build data using supplied 'character' probe vector
            if(class(probes) == "character"){

              data <- t(object@exprs[probes, ])
            }

            # Set up and perform eBayes
            design <- as.matrix(ifelse(object@annot$defineCase == "Case", 1, 0))
            colnames(design) <- "CaseVCont"
            fit <- limma::lmFit(object@exprs[probes, ], design)
            ebaye <- limma::ebayes(fit, ...)

            # Sort probes
            vals <- data.frame("probe" = rownames(ebaye$p.value), "p.value" = ebaye$p.value[, 1])
            final <- as.character(vals[order(vals$p.value), "probe"])

            array <- new("ExprsBinary",
                         exprs = object@exprs[final,],
                         annot = object@annot,
                         preFilter = append(object@preFilter, list(final)),
                         reductionModel = append(object@reductionModel, list(NA))
            )

            return(array)
          }
)

#' @rdname fs
#' @section Methods (by generic):
#' \code{fsMrme:} Method to perform mRMR feature selection using mRMRe::mRMR.classic.
#' @importFrom mRMRe mRMR.data mRMR.classic solutions featureNames
#' @export
setMethod("fsMrmre", "ExprsBinary",
          function(object, probes, ...){ # args to mRMR.classic

            args <- as.list(substitute(list(...)))[-1]

            if(!"target_indices" %in% names(args)){

              cat("Setting 'target_indices' to 1 (default behavior, override explicitly)...\n")
              args <- append(args, list("target_indices" = 1))
            }

            if(!"feature_count" %in% names(args)){

              cat("Setting 'feature_count' to 64 (default behavior, override explicitly)...\n")
              args <- append(args, list("feature_count" = 64))
            }

            # Convert 'numeric' probe argument to 'character' probe vector
            if(class(probes) == "numeric"){

              if(probes == 0) probes <- nrow(object@exprs)
              probes <- rownames(object@exprs[1:probes, ])
              data <- t(object@exprs[probes, ])
            }

            # Build data using supplied 'character' probe vector
            if(class(probes) == "character"){

              data <- t(object@exprs[probes, ])
            }

            # Set up "make.names" key for improper @exprs row.names
            key <- data.frame("old" = colnames(data), "new" = make.names(colnames(data)))

            # Set up and perform mRMR
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
            final <- merge(data.frame("new" = final), key)$old

            array <- new("ExprsBinary",
                         exprs = object@exprs[final,],
                         annot = object@annot,
                         preFilter = append(object@preFilter, list(final)),
                         reductionModel = append(object@reductionModel, list(NA))
            )

            return(array)

          }
)
