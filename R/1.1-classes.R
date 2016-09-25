###########################################################
### ExprsArray class

#' An S4 class to store feature and annotation data
#'
#' @slot exprs A matrix. Stores the feature data.
#' @slot annot A data.frame. Stores the annotation data.
#' @slot preFilter Typically a list. Stores feature selection history.
#' @slot reductionModel Typically a list. Stores dimension reduction history.
#'
#' @seealso
#' \code{\link{ExprsArray-class}}\cr
#' \code{\link{ExprsModel-class}}\cr
#' \code{\link{ExprsPipeline-class}}\cr
#' \code{\link{ExprsEnsemble-class}}\cr
#' \code{\link{ExprsPredict-class}}
#' @export
setClass("ExprsArray",
         slots = c(
           exprs = "matrix",
           annot = "data.frame",
           preFilter = "ANY",
           reductionModel = "ANY"
         )
)

#' An S4 class to store feature and annotation data
#'
#' An \code{ExprsArray} sub-class for data with binary class labels.
#'
#' @export
setClass("ExprsBinary",
         contains = "ExprsArray"
)

#' An S4 class to store feature and annotation data
#'
#' An \code{ExprsArray} sub-class for data with many class labels.
#'
#' @export
setClass("ExprsMulti",
         contains = "ExprsArray"
)

###########################################################
### ExprsModel class

#' An S4 class to store the classification model
#'
#' @slot preFilter Typically a list. Stores feature selection history.
#' @slot reductionModel Typically a list. Stores dimension reduction history.
#' @slot mach Typically an S4 class. Stores the classification model.
#'
#' @seealso
#' \code{\link{ExprsArray-class}}\cr
#' \code{\link{ExprsModel-class}}\cr
#' \code{\link{ExprsPipeline-class}}\cr
#' \code{\link{ExprsEnsemble-class}}\cr
#' \code{\link{ExprsPredict-class}}
#' @export
setClass("ExprsModel",
         slots = c(
           preFilter = "ANY",
           reductionModel = "ANY",
           mach = "ANY"
         )
)

#' An S4 class to store the classification model
#'
#' An \code{ExprsModel} sub-class for dichotomous classifiers.
#'
#' @export
setClass("ExprsMachine",
         contains = "ExprsModel"
)

#' An S4 class to store the classification model
#'
#' An \code{ExprsModel} sub-class for multi-class classifiers.
#'
#' @export
setClass("ExprsModule",
         contains = "ExprsModel"
)

###########################################################
### ExprsPipeline class

#' An S4 class to store models built during high-throughput learning
#'
#' @slot summary Typically a data.frame. Stores the parameters and
#'  performances for classification models.
#' @slot machs Typically a list. Stores the classification models
#'  referenced in \code{summary} slot.
#'
#' @seealso
#' \code{\link{ExprsArray-class}}\cr
#' \code{\link{ExprsModel-class}}\cr
#' \code{\link{ExprsPipeline-class}}\cr
#' \code{\link{ExprsEnsemble-class}}\cr
#' \code{\link{ExprsPredict-class}}
#' @export
setClass("ExprsPipeline",
         slots = c(
           summary = "ANY",
           machs = "ANY")
)

###########################################################
### ExprsEnsemble class

#' An S4 class to store multiple classification models
#'
#' @slot machs Typically a list. Stores the classification models
#'  referenced in \code{summary} slot.
#'
#' @seealso
#' \code{\link{ExprsArray-class}}\cr
#' \code{\link{ExprsModel-class}}\cr
#' \code{\link{ExprsPipeline-class}}\cr
#' \code{\link{ExprsEnsemble-class}}\cr
#' \code{\link{ExprsPredict-class}}
#' @export
setClass("ExprsEnsemble",
         slots = c(
           machs = "ANY"
         )
)

###########################################################
### ExprsPredict class

#' An S4 class to store class predictions
#'
#' @slot pred A factor. Stores class predictions as an unambiguous
#'  class assignment.
#' @slot decision.values Typically a matrix. Stores class predictions
#'  as a decision value.
#' @slot probability Typically a matrix. Stores class predictions
#'  as a probability.
#' @slot actual Typically a factor. Stores known class labels.
#'  Used by \code{\link{calcStats}}.
#'
#' @seealso
#' \code{\link{ExprsArray-class}}\cr
#' \code{\link{ExprsModel-class}}\cr
#' \code{\link{ExprsPipeline-class}}\cr
#' \code{\link{ExprsEnsemble-class}}\cr
#' \code{\link{ExprsPredict-class}}
#' @export
setClass("ExprsPredict",
         slots = c(
           pred = "factor",
           decision.values = "ANY",
           probability = "ANY",
           actual = "ANY"
         )
)
