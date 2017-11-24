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
#' \code{\link{ExprsPredict-class}}\cr
#' \code{\link{RegrsPredict-class}}
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
#' An \code{ExprsArray} sub-class for data with binary class outcomes.
#'
#' @export
setClass("ExprsBinary",
         contains = "ExprsArray"
)

#' An S4 class to store feature and annotation data
#'
#' An \code{ExprsArray} sub-class for data with multiple class outcomes.
#'
#' @export
setClass("ExprsMulti",
         contains = "ExprsArray"
)

#' An S4 class to store feature and annotation data
#'
#' An \code{ExprsArray} sub-class for data with continuous outcomes.
#'
#' @export
setClass("RegrsArray",
         contains = "ExprsArray"
)

#' An S4 class to store the model
#'
#' @slot preFilter Typically a list. Stores feature selection history.
#' @slot reductionModel Typically a list. Stores dimension reduction history.
#' @slot mach Typically an S4 class. Stores the model.
#'
#' @seealso
#' \code{\link{ExprsArray-class}}\cr
#' \code{\link{ExprsModel-class}}\cr
#' \code{\link{ExprsPipeline-class}}\cr
#' \code{\link{ExprsEnsemble-class}}\cr
#' \code{\link{ExprsPredict-class}}\cr
#' \code{\link{RegrsPredict-class}}
#' @export
setClass("ExprsModel",
         slots = c(
           preFilter = "ANY",
           reductionModel = "ANY",
           mach = "ANY"
         )
)

#' An S4 class to store the model
#'
#' An \code{ExprsModel} sub-class for dichotomous models.
#'
#' @export
setClass("ExprsMachine",
         contains = "ExprsModel"
)

#' An S4 class to store the model
#'
#' An \code{ExprsModel} sub-class for multi-class models.
#'
#' @export
setClass("ExprsModule",
         contains = "ExprsModel"
)

#' An S4 class to store the model
#'
#' An \code{ExprsModel} sub-class for continuous outcome models.
#'
#' @export
setClass("RegrsModel",
         contains = "ExprsModel"
)

#' An S4 class to store models built during high-throughput learning
#'
#' @slot summary Typically a data.frame. Stores the parameters and
#'  performances for the models.
#' @slot machs Typically a list. Stores the models
#'  referenced in \code{summary} slot.
#'
#' @seealso
#' \code{\link{ExprsArray-class}}\cr
#' \code{\link{ExprsModel-class}}\cr
#' \code{\link{ExprsPipeline-class}}\cr
#' \code{\link{ExprsEnsemble-class}}\cr
#' \code{\link{ExprsPredict-class}}\cr
#' \code{\link{RegrsPredict-class}}
#' @export
setClass("ExprsPipeline",
         slots = c(
           summary = "ANY",
           machs = "ANY")
)

#' An S4 class to store multiple models
#'
#' @slot machs Typically a list. Stores the models.
#'
#' @seealso
#' \code{\link{ExprsArray-class}}\cr
#' \code{\link{ExprsModel-class}}\cr
#' \code{\link{ExprsPipeline-class}}\cr
#' \code{\link{ExprsEnsemble-class}}\cr
#' \code{\link{ExprsPredict-class}}\cr
#' \code{\link{RegrsPredict-class}}
#' @export
setClass("ExprsEnsemble",
         slots = c(
           machs = "ANY"
         )
)

#' An S4 class to store model predictions
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
#' \code{\link{ExprsPredict-class}}\cr
#' \code{\link{RegrsPredict-class}}
#' @export
setClass("ExprsPredict",
         slots = c(
           pred = "factor",
           decision.values = "ANY",
           probability = "ANY",
           actual = "ANY"
         )
)

#' An S4 class to store model predictions
#'
#' @slot pred Any. Stores predicted outcome.
#' @slot actual Any. Stores actual outcome.
#'  Used by \code{\link{calcStats}}.
#'
#' @seealso
#' \code{\link{ExprsArray-class}}\cr
#' \code{\link{ExprsModel-class}}\cr
#' \code{\link{ExprsPipeline-class}}\cr
#' \code{\link{ExprsEnsemble-class}}\cr
#' \code{\link{ExprsPredict-class}}\cr
#' \code{\link{RegrsPredict-class}}
#' @export
setClass("RegrsPredict",
         slots = c(
           pred = "ANY",
           actual = "ANY"
         )
)
