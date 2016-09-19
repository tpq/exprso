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
#'  parameter space in a high-throughput manner, see \code{pl} methods.
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
#' @inheritParams fs
#' @param object Specifies the \code{ExprsArray} object to use as a training set
#'  for classification.
#' @return Returns an \code{ExprsModel} object.
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
           function(object, probes = 0, ...) standardGeneric("buildNB")
)

#' @rdname build
#' @export
setGeneric("buildLDA",
           function(object, probes = 0, ...) standardGeneric("buildLDA")
)

#' @rdname build
#' @export
setGeneric("buildSVM",
           function(object, probes = 0, ...) standardGeneric("buildSVM")
)

#' @rdname build
#' @export
setGeneric("buildANN",
           function(object, probes = 0, ...) standardGeneric("buildANN")
)

#' @rdname build
#' @export
setGeneric("buildRF",
           function(object, probes = 0, ...) standardGeneric("buildRF")
)

#' @rdname build
#' @export
setGeneric("buildDNN",
           function(object, probes = 0, ...) standardGeneric("buildDNN")
)

###########################################################
### Build classifier

#' Workhorse for build Methods
#'
#' Used as a back-end wrapper for creating new build methods.
#'
#' @inheritParams build
#' @param uniqueFx A function call unique to that fs method.
#' @return Returns an \code{ExprsModel} object.
#'
#' @export
build. <- function(object, probes, uniqueFx, ...){

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
  labels <- factor(object@annot[rownames(data), "defineCase"], levels = c("Control", "Case"))
  model <- do.call("uniqueFx", list(data, labels, ...))

  # Carry through and append fs history as stored in the ExprsArray object
  # NOTE: length(ExprsMachine@preFilter) > length(ExprsArray@preFilter)
  machine <- new("ExprsMachine",
                 preFilter = append(object@preFilter, list(probes)),
                 reductionModel = append(object@reductionModel, list(NA)),
                 mach = model
  )
}

#' @rdname build
#' @section Methods (by generic):
#' \code{buildNB:} Method to build classifiers using e1071::naiveBayes.
#' @importFrom e1071 naiveBayes
#' @export
setMethod("buildNB", "ExprsBinary",
          function(object, probes, ...){ # args to naiveBayes

            build.(object, probes,
                   uniqueFx = function(data, labels, ...){

                     # Perform naiveBayes via ~ method
                     args <- getArgs(...)
                     df <- data.frame(data, "defineCase" = labels)
                     args <- append(list("formula" = defineCase ~ ., "data" = df), args)
                     do.call(e1071::naiveBayes, args)
                   }, ...)
          }
)

#' @rdname build
#' @section Methods (by generic):
#' \code{buildLDA:} Method to build classifiers using MASS::lda.
#' @importFrom MASS lda
#' @export
setMethod("buildLDA", "ExprsBinary",
          function(object, probes, ...){ # args to lda

            build.(object, probes,
                   uniqueFx = function(data, labels, ...){

                     # Perform linear discriminant analysis via ~ method
                     args <- getArgs(...)
                     df <- data.frame(data, "defineCase" = labels)
                     args <- append(list("formula" = defineCase ~ ., "data" = df), args)
                     do.call(MASS::lda, args)
                   }, ...)
          }
)

#' @rdname build
#' @section Methods (by generic):
#' \code{buildSVM:} Method to build classifiers using e1071::svm.
#' @importFrom e1071 svm
#' @export
setMethod("buildSVM", "ExprsBinary",
          function(object, probes, ...){ # args to svm

            build.(object, probes,
                   uniqueFx = function(data, labels, ...){

                     # Perform SVM via ~ method (permits plotting)
                     args <- getArgs(...)
                     args <- forceArg("probability", TRUE, args)
                     args <- forceArg("cross", 0, args)
                     df <- data.frame(data, "defineCase" = labels)
                     args <- append(list("formula" = defineCase ~ ., "data" = df), args)
                     do.call(e1071::svm, args)
                   }, ...)
          }
)

#' @rdname build
#' @section Methods (by generic):
#' \code{buildANN:} Method to build classifiers using nnet::nnet.
#' @importFrom nnet nnet
#' @export
setMethod("buildANN", "ExprsBinary",
          function(object, probes, ...){ # args to nnet

            build.(object, probes,
                   uniqueFx = function(data, labels, ...){

                     # Perform ANN via ~ method
                     args <- getArgs(...)
                     args <- defaultArg("size", 1, args)
                     args <- defaultArg("range", 1/max(abs(as.vector(data))), args)
                     args <- defaultArg("decay", 0.5, args)
                     args <- defaultArg("maxit", 1000, args)
                     df <- data.frame(data, "defineCase" = labels)
                     args <- append(list("formula" = defineCase ~ ., "data" = df), args)
                     do.call(nnet::nnet, args)
                   }, ...)
          }
)

#' @rdname build
#' @section Methods (by generic):
#' \code{buildRF:} Method to build classifiers using randomForest::randomForest.
#' @importFrom randomForest randomForest
#' @export
setMethod("buildRF", "ExprsBinary",
          function(object, probes, ...){ # args to randomForest

            build.(object, probes,
                   uniqueFx = function(data, labels, ...){

                     # Perform RF via ~ method
                     args <- getArgs(...)
                     df <- data.frame(data, "defineCase" = labels)
                     args <- append(list("formula" = defineCase ~ ., "data" = df), args)
                     do.call(randomForest::randomForest, args)
                   }, ...)
          }
)

#' @rdname build
#' @section Methods (by generic):
#' \code{buildDNN:} Method to build feed-forward networks using h2o::h2o.deeplearning.
#' @importFrom utils write.csv
#' @importFrom h2o h2o.init h2o.importFile h2o.deeplearning
#' @export
setMethod("buildDNN", "ExprsBinary",
          function(object, probes, ...){ # args to h2o.deeplearning

            args <- getArgs(...)

            build.(object, probes,
                   uniqueFx = function(data, labels, ...){

                     # Initialize h2o (required)
                     localH20 <- h2o::h2o.init()

                     # Rename data based on order supplied by probes
                     colnames(data) <- paste0("id", 1:ncol(data))
                     labels <- object@annot[rownames(data), "defineCase"]
                     df <- data.frame(data, "defineCase" = labels)

                     # Import data as H2OFrame via a temporary csv
                     tempFile <- tempfile(fileext = ".csv")
                     utils::write.csv(df, tempFile)
                     h2o.data <- h2o::h2o.importFile(path = tempFile, destination_frame = "h2o.data")

                     # Prepare arguments and build classifier
                     args <- getArgs(...)
                     args <- append(list("x" = colnames(data), "y" = "defineCase",
                                         "training_frame" = h2o.data), args)
                     do.call(h2o::h2o.deeplearning, args)
                   }, ...)
          }
)
