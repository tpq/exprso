#' Build Ensemble
#'
#' \code{buildEnsemble} builds an ensemble from \code{ExprsModel} or
#'  \code{ExprsPipeline} objects. See Details.
#'
#' This function can combine any number of model objects into an ensemble.
#'  These models do not necessarily have to derive from the same \code{build}
#'  method. In this way, it works like \code{\link{conjoin}} method.
#'
#' This function can also build an ensemble from pipeline objects. It does
#'  this by calling \code{\link{pipeFilter}}, then aggregating those results
#'  into an ensemble. As an adjunct to this method, consider first combining
#'  multiple pipeline objects together with \code{\link{conjoin}}.
#'
#' @param object An \code{ExprsModel} or \code{ExprsPipeline} object.
#' @param ... Additional \code{ExprsModel} objects to use in the ensemble.
#'  Argument applies to the \code{\link{ExprsModel-class}} method only.
#' @inheritParams pipeFilter
#' @return An \code{\link{ExprsEnsemble-class}} object.
#' @export
setGeneric("buildEnsemble",
           function(object, ...) standardGeneric("buildEnsemble")
)

#' @describeIn buildEnsemble Method to build ensemble from \code{ExprsModel} objects.
#' @export
setMethod("buildEnsemble", "ExprsModel",
          function(object, ...){ # args to include additional ExprsMachine objects

            conjoin(object, ...)
          }
)

#' @describeIn buildEnsemble Method to build ensemble from \code{ExprsPipeline} objects.
#' @export
setMethod("buildEnsemble", "ExprsPipeline",
          function(object, colBy = 0, how = 0, gate = 0, top = 0){

            if(!identical(colBy, 0)){

              object <- pipeFilter(object, colBy = colBy, how = how, gate = gate, top = top)
            }

            new("ExprsEnsemble",
                machs = unlist(object@machs)
            )
          }
)
