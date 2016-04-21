###########################################################
### Define generic functions

#' @name build
#' @rdname build
#'
#' @title Build Classifiers
#'
#' @description A collection of functions to build classifiers.
#'
#' @details
#'
#' These \code{build} methods construct a single classifier given an \code{ExprsArray}
#'  object and a set of parameters. This function returns an \code{ExprsModel} object.
#'  In the case of binary classification, these methods use an \code{ExprsBinary}
#'  object and return an \code{ExprsMachine} object. In the case of multi-class
#'  classification, these methods use an \code{ExprsMulti} object and return an
#'  \code{ExprsModule} object. In the case of multi-class classification, these methods
#'  harness the \code{\link{doMulti}} function to perform "1 vs. all" classifier
#'  construction. In the setting of four class labels, a single \code{build} call
#'  will return four classifiers that work in concert to make a single prediction
#'  of an unlabelled subject. For building multiple classifiers across a vast
#'  parameter space in a high-throughput manner, see \code{\link{pl}} methods.
#'
#' Like \code{\link{fs}} methods, \code{build} methods have a \code{probes} argument
#'  which allows the user to specify which features to feed INTO the classifier
#'  build. This effectively provides the user with one last opportunity to subset
#'  the feature space based on prior feature selection or dimension reduction.
#'  For all build methods, \code{@@preFilter} and \code{@@reductionModel} will
#'  get passed along to the resultant \code{ExprsModel} object, again ensuring
#'  that any test or validation sets will undergo the same feature selection and
#'  dimension reduction in the appropriate steps when deploying the classifier.
#'  Set \code{probes = 0} to pass all features through a \code{build} method.
#'  See \code{\link{modHistory}} to read more about feature selection history.
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

#' @rdname build
#' @export
setGeneric("buildNB",
           function(object, ...) standardGeneric("buildNB")
)

#' @rdname build
#' @export
setGeneric("buildLDA",
           function(object, ...) standardGeneric("buildLDA")
)

#' @rdname build
#' @export
setGeneric("buildSVM",
           function(object, ...) standardGeneric("buildSVM")
)

#' @rdname build
#' @export
setGeneric("buildANN",
           function(object, ...) standardGeneric("buildANN")
)

#' @rdname build
#' @export
setGeneric("buildRF",
           function(object, ...) standardGeneric("buildRF")
)

###########################################################
### Build classifier

#' @rdname build
#' @section Methods (by generic):
#' \code{buildNB:} Method to build classifiers using e1071::naiveBayes.
#'
#' @inheritParams fs
#' @return Returns an \code{ExprsModel} object.
#'
#' @importFrom e1071 naiveBayes
#' @export
setMethod("buildNB", "ExprsBinary",
          function(object, probes, ...){ # args to naiveBayes

            args <- as.list(substitute(list(...)))[-1]

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

            # Perform naiveBayes via ~ method
            df <- data.frame(data, "defineCase" = labels)
            args <- append(list("formula" = defineCase ~ ., "data" = df), args)
            model <- do.call(e1071::naiveBayes, args)

            # Carry through and append fs history as stored in the ExprsArray object
            # NOTE: length(ExprsMachine@preFilter) > length(ExprsArray@preFilter)
            machine <- new("ExprsMachine",
                           preFilter = append(object@preFilter, list(probes)),
                           reductionModel = append(object@reductionModel, list(NA)),
                           mach = model
            )

            return(machine)
          }
)

#' @rdname build
#' @section Methods (by generic):
#' \code{buildLDA:} Method to build classifiers using MASS::lda.
#' @importFrom MASS lda
#' @export
setMethod("buildLDA", "ExprsBinary",
          function(object, probes, ...){ # args to lda

            args <- as.list(substitute(list(...)))[-1]

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

            # Perform linear discriminant analysis via ~ method
            df <- data.frame(data, "defineCase" = labels)
            args <- append(list("formula" = defineCase ~ ., "data" = df), args)
            model <- do.call(MASS::lda, args)

            # Carry through and append fs history as stored in the ExprsArray object
            # NOTE: length(ExprsMachine@preFilter) > length(ExprsArray@preFilter)
            machine <- new("ExprsMachine",
                           preFilter = append(object@preFilter, list(probes)),
                           reductionModel = append(object@reductionModel, list(NA)),
                           mach = model
            )

            return(machine)
          }
)

#' @rdname build
#' @section Methods (by generic):
#' \code{buildSVM:} Method to build classifiers using e1071::svm.
#' @importFrom e1071 svm
#' @export
setMethod("buildSVM", "ExprsBinary",
          function(object, probes, ...){ # args to svm

            args <- as.list(substitute(list(...)))[-1]

            if(!"probability" %in% names(args)){

              args <- append(args, list("probability" = TRUE))
            }

            # Mandate probability = TRUE
            if(args$probability == FALSE){

              cat("Uh oh! This function requires 'probability' = TRUE. Setting 'probability' to TRUE...\n")
              args$probability <- TRUE
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

            # Set labels as factor
            labels <- factor(object@annot[rownames(data), "defineCase"], levels = c("Control", "Case"))

            # Perform SVM via ~ method (permits plotting)
            df <- data.frame(data, "defineCase" = labels)
            args <- append(list("formula" = defineCase ~ ., "data" = df), args)
            model <- do.call(e1071::svm, args)

            # Carry through and append fs history as stored in the ExprsArray object
            # NOTE: length(ExprsMachine@preFilter) > length(ExprsArray@preFilter)
            machine <- new("ExprsMachine",
                           preFilter = append(object@preFilter, list(probes)),
                           reductionModel = append(object@reductionModel, list(NA)),
                           mach = model
            )

            return(machine)
          }
)

#' @rdname build
#' @section Methods (by generic):
#' \code{buildANN:} Method to build classifiers using nnet::nnet.
#' @importFrom nnet nnet
#' @export
setMethod("buildANN", "ExprsBinary",
          function(object, probes, ...){

            args <- as.list(substitute(list(...)))[-1]

            if(!"size" %in% names(args)){

              cat("Setting 'size' to 1 (default behavior, override explicitly)...\n")
              args <- append(args, list("size" = 1))
            }

            if(!"range" %in% names(args)){

              cat("Setting 'range' to 1/max(|x|) (default behavior, override explicitly)...\n")
              args <- append(args, list("range" = 1/max(abs(as.vector(object@exprs)))))
            }

            if(!"decay" %in% names(args)){

              cat("Setting 'decay' to 0.5 (default behavior, override explicitly)...\n")
              args <- append(args, list("decay" = 0.5))
            }

            if(!"maxit" %in% names(args)){

              cat("Setting 'maxit' to 1000 (default behavior, override explicitly)...\n")
              args <- append(args, list("maxit" = 1000))
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

            # Set labels as factor
            labels <- factor(object@annot[rownames(data), "defineCase"], levels = c("Control", "Case"))

            # Perform ANN via ~ method
            df <- data.frame(data, "defineCase" = labels)
            args <- append(list("formula" = defineCase ~ ., "data" = df), args)
            model <- do.call(nnet::nnet, args)

            # Carry through and append fs history as stored in the ExprsArray object
            # NOTE: length(ExprsMachine@preFilter) > length(ExprsArray@preFilter)
            machine <- new("ExprsMachine",
                           preFilter = append(object@preFilter, list(probes)),
                           reductionModel = append(object@reductionModel, list(NA)),
                           mach = model
            )

            return(machine)
          }
)

#' @rdname build
#' @section Methods (by generic):
#' \code{buildRF:} Method to build classifiers using randomForest::randomForest.
#' @importFrom randomForest randomForest
#' @export
setMethod("buildRF", "ExprsBinary",
          function(object, probes, ...){ # args to randomForest

            args <- as.list(substitute(list(...)))[-1]

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

            # Perform RF via ~ method
            df <- data.frame(data, "defineCase" = labels)
            args <- append(list("formula" = defineCase ~ ., "data" = df), args)
            model <- do.call(randomForest::randomForest, args)

            # Carry through and append fs history as stored in the ExprsArray object
            # NOTE: length(ExprsMachine@preFilter) > length(ExprsArray@preFilter)
            machine <- new("ExprsMachine",
                           preFilter = append(object@preFilter, list(probes)),
                           reductionModel = append(object@reductionModel, list(NA)),
                           mach = model
            )

            return(machine)
          }
)
