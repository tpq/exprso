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
#'  to prioritize which features to include during classifier construction. Although there exists
#'  many feature selection methods, this package provides wrappers for some of the most popular ones.
#'  Each wrapper (1) pre-processes the \code{ExprsArray} input, (2) performs the feature selection,
#'  and (3) returns an \code{ExprsArray} output with an updated feature selection history.
#'  You can use, in tandem, any number of feature selection methods, and in any order.
#'
#' For all feature selection methods, \code{@@preFilter} and \code{@@reductionModel} stores the
#'  feature selection and dimension reduction history, respectively. This history gets passed
#'  along to prepare the test or validation set during model deployment, ensuring that these
#'  sets undergo the same feature selection and dimension reduction as the training set.
#'
#' Under the scenarios where users plan to apply multiple feature selection or dimension
#'  reduction steps, the \code{top} argument manages which features (e.g., gene expression values)
#'  to send through each feature selection or dimension reduction procedure. For \code{top},
#'  a numeric scalar indicates the number of top features to use, while a character vector
#'  indicates specifically which features to use. In this way, the user sets which features
#'  to feed INTO the \code{fs} method (NOT which features the user expects OUT). The example
#'  below shows how to apply dimension reduction to the top 50 features as selected by the
#'  Student's t-test. Set \code{top = 0} to pass all features through an \code{fs} method.
#'
#' Note that not all feature selection methods will generalize to multi-class data.
#'  A feature selection method will fail when applied to an \code{ExprsMulti} object
#'  unless that feature selection method has an \code{ExprsMulti} method.
#'
#' Note that \code{fsMrmre} crashes when supplied a very large \code{feature_count} argument
#'  owing to its implementation in the imported package \code{mRMRe}.
#'
#' @param object Specifies the \code{ExprsArray} object to undergo feature selection.
#' @param top A numeric scalar or character vector. A numeric scalar indicates
#'  the number of top features that should undergo feature selection. A character vector
#'  indicates specifically which features by name should undergo feature selection.
#'  Set \code{top = 0} to include all features. A numeric vector can also be used
#'  to indicate specific features by location, similar to a character vector.
#' @param how A character string. Toggles between the sub-routines "t.test" and
#'  "ks.test". Argument applies to \code{fsStats} only.
#' @param ... Arguments passed to the respective wrapped function.
#'
#' @return Returns an \code{ExprsArray} object.
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
#' library(golubEsets)
#' data(Golub_Merge)
#' array <- arrayEset(Golub_Merge, colBy = "ALL.AML", include = list("ALL", "AML"))
#' array <- modFilter(array, 20, 16000, 500, 5) # pre-filter Golub ala Deb 2003
#' array <- modTransform(array) # lg transform
#' array <- modNormalize(array, c(1, 2)) # normalize gene and subject vectors
#' arrays <- splitSample(array, percent.include = 67)
#' array.train <- fsStats(arrays[[1]], top = 0, how = "t.test")
#' array.train <- fsPrcomp(array.train, top = 50)
#' mach <- buildSVM(array.train, top = 5, kernel = "linear", cost = 1)
#' }
NULL

#' @rdname fs
#' @export
setGeneric("fsSample",
           function(object, top = 0, ...) standardGeneric("fsSample")
)

#' @rdname fs
#' @export
setGeneric("fsNULL",
           function(object, top = 0, ...) standardGeneric("fsNULL")
)

#' @rdname fs
#' @export
setGeneric("fsANOVA",
           function(object, top = 0, ...) standardGeneric("fsANOVA")
)

#' @rdname fs
#' @export
setGeneric("fsStats",
           function(object, top = 0, ...) standardGeneric("fsStats")
)

#' @rdname fs
#' @export
setGeneric("fsPrcomp",
           function(object, top = 0, ...) standardGeneric("fsPrcomp")
)

#' @rdname fs
#' @export
setGeneric("fsPathClassRFE",
           function(object, top = 0, ...) standardGeneric("fsPathClassRFE")
)

#' @rdname fs
#' @export
setGeneric("fsEbayes",
           function(object, top = 0, ...) standardGeneric("fsEbayes")
)

#' @rdname fs
#' @export
setGeneric("fsMrmre",
           function(object, top = 0, ...) standardGeneric("fsMrmre")
)

###########################################################
### Select features

#' Workhorse for fs Methods
#'
#' Used as a back-end wrapper for creating new fs methods.
#'
#' If the uniqueFx returns a character vector, it is assumed
#'  that the fs method is for feature selection only. If the
#'  uniqueFx returns a list, it is assumed that the fs method
#'  is a reduction model method only.
#'
#' @inheritParams fs
#' @param uniqueFx A function call unique to that fs method.
#' @return Returns an \code{ExprsArray} object.
#'
#' @export
fs. <- function(object, top, uniqueFx, ...){

  if(class(top) == "numeric"){

    if(length(top) == 1){

      if(top > nrow(object@exprs)) top <- 0
      if(top == 0) top <- nrow(object@exprs)
      top <- rownames(object@exprs[1:top, ])

    }else{

      top <- rownames(object@exprs[top, ])
    }
  }

  data <- t(object@exprs[top, ])
  final <- do.call("uniqueFx", list(data, top, ...))

  if(class(final) == "character"){

    array <- new(class(object), exprs = object@exprs[final,], annot = object@annot,
                 preFilter = append(object@preFilter, list(final)),
                 reductionModel = append(object@reductionModel, list(NA))
    )

  }else if(class(final) == "list"){

    array <- new(class(object), exprs = final[[1]], annot = object@annot,
                 preFilter = append(object@preFilter, list(top)),
                 reductionModel = append(object@reductionModel, list(final[[2]]))
    )

  }else{

    stop("Uh oh! DEBUG ERROR: 002")
  }

  return(array)
}

#' @rdname fs
#' @section Methods (by generic):
#' \code{fsSample:} Method to perform random feature selection using base::sample.
#' @export
setMethod("fsSample", "ExprsArray",
          function(object, top, ...){ # args to sample

            fs.(object, top,
                uniqueFx = function(data, top, ...){

                  sample(top, ...)
                }, ...)
          }
)

#' @rdname fs
#' @section Methods (by generic):
#' \code{fsNULL:} Method to perform a NULL feature selection and return input unaltered.
#' @export
setMethod("fsNULL", "ExprsArray",
          function(object, top, ...){ # args to NULL

            fs.(object, top,
                uniqueFx = function(data, top, ...){

                  top
                }, ...)
          }
)

#' @rdname fs
#' @section Methods (by generic):
#' \code{fsANOVA:} Method to perform ANOVA feature selection using stats::aov.
#' @importFrom stats aov
#' @export
setMethod("fsANOVA", "ExprsArray",
          function(object, top, ...){ # args to aov

            fs.(object, top,
                uniqueFx = function(data, top, ...){

                  # Perform an ANOVA for each feature in data
                  df <- data.frame(data, "label" = object@annot[rownames(data), "defineCase"])
                  p <- vector("numeric", length(top))
                  for(i in 1:length(top)){

                    formula <- stats::as.formula(paste(top[i], "~", "label"))
                    fit <- stats::aov(formula, data = df, ...)
                    p[i] <- summary(fit)[[1]][1, "Pr(>F)"]
                  }

                  top[order(p)]
                }, ...)
          }
)

#' @rdname fs
#' @section Methods (by generic):
#' \code{fsStats:} Method to perform statistics based feature selection using stats::t.test and others.
#' @importFrom stats t.test ks.test
#' @export
setMethod("fsStats", "ExprsBinary",
          function(object, top = 0, how = c("t.test", "ks.test"), ...){ # args to t.test or ks.test

            # Choose first "how" from default argument vector
            how <- how[1]

            if(how == "t.test"){

              fs.(object, top,
                  uniqueFx = function(data, top, ...){

                    # Prepare data for statistical tests
                    cases <- object@annot$defineCase %in% "Case"
                    conts <- object@annot$defineCase %in% "Control"
                    p <- vector("numeric", length(top))

                    for(i in 1:length(top)){

                      tryCatch(
                        {
                          p[i] <- t.test(object@exprs[top[i], cases],
                                         object@exprs[top[i], conts], ...)$p.value

                        }, error = function(e){

                          cat("fsStats failed for feature: ", top[i], ". Setting p(x)=1...\n")
                          p[i] <- 1
                        })
                    }

                    top[order(p)]
                  }, ...)

            }else if(how == "ks.test"){

              fs.(object, top,
                  uniqueFx = function(data, top, ...){

                    # Prepare data for statistical tests
                    cases <- object@annot$defineCase %in% "Case"
                    conts <- object@annot$defineCase %in% "Control"
                    p <- vector("numeric", length(top))

                    for(i in 1:length(top)){

                      tryCatch(
                        {
                          p[i] <- ks.test(object@exprs[top[i], cases],
                                          object@exprs[top[i], conts], ...)$p.value

                        }, error = function(e){

                          cat("fsStats failed for feature: ", top[i], ". Setting p(x)=1...\n")
                          p[i] <- 1
                        })
                    }

                    top[order(p)]
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
          function(object, top, ...){ # args to prcomp

            fs.(object, top,
                uniqueFx = function(data, top, ...){

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
          function(object, top, ...){ # args to fit.rfe

            fs.(object, top,
                uniqueFx = function(data, top, ...){

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
#' @export
setMethod("fsEbayes", "ExprsBinary",
          function(object, top, ...){ # args to ebayes

            if(!requireNamespace("limma", quietly = TRUE)){
              stop("limma needed for this function to work. Please install it.",
                   call. = FALSE)
            }

            fs.(object, top,
                uniqueFx = function(data, top, ...){

                  design <- as.matrix(ifelse(object@annot$defineCase == "Case", 1, 0))
                  colnames(design) <- "CaseVCont"
                  fit <- limma::lmFit(object@exprs[top, ], design)
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
          function(object, top, ...){ # args to mRMR.classic

            fs.(object, top,
                uniqueFx = function(data, top, ...){

                  args <- getArgs(...)
                  args <- defaultArg("target_indices", 1, args)
                  args <- defaultArg("feature_count", 64, args)

                  # Set up "make.names" key for improper @exprs row.names
                  key <- data.frame("old" = colnames(data), "new" = make.names(colnames(data)))

                  labels <- as.numeric(object@annot$defineCase == "Case")
                  mRMRdata <- mRMRe::mRMR.data(data = data.frame(labels, data))
                  args <- append(list("data" = mRMRdata), args)
                  mRMRout <- do.call(mRMRe::mRMR.classic, args)

                  # Sort features
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
