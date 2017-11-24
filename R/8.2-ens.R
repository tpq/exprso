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

#' @rdname exprso-predict
#' @export
setMethod("predict", "ExprsEnsemble",
          function(object, array, how = "probability", verbose = TRUE){

            classCheck(array, "ExprsArray",
                       "You can only deploy a model on the results of ?exprso.")

            # Deploy each machine in @machs on the provided ExprsArray
            results <- lapply(object@machs, function(mach) predict(mach, array, verbose = verbose))

            if(class(array) == "ExprsMulti" |
               (class(array) == "ExprsBinary" & how == "majority")){

              # Cast a vote (1 for Case, -1 for Control)
              votes <- lapply(results, function(result) ifelse(result@pred == "Case", 1, -1))
              final <- rowSums(data.frame(votes, row.names = names(results[[1]]@pred)))

              # Randomly assign ties as "Case" or "Control" based on the proportion of cases
              case <- sum(array@annot$defineCase == "Case") / nrow(array@annot)
              final[final == 0] <- sample(c(1, -1),
                                          size = length(final[final == 0]),
                                          replace = TRUE,
                                          prob = c(case, 1 - case))

              # Clean up pred
              pred <- factor(ifelse(final > 0, "Case", "Control"), levels = c("Control", "Case"))
              px <- NULL
              dv <- NULL

              final <- new("ExprsPredict",
                           pred = pred, decision.values = dv, probability = px,
                           actual = array@annot$defineCase)

            }else if(class(array) == "ExprsBinary" & how == "probability"){

              # Calculate average probabilities
              pxs <- lapply(results, function(result) result@probability)
              Case <- lapply(pxs, function(px) px[, "Case"])
              Case <- rowMeans(as.data.frame(Case))
              Control <- lapply(pxs, function(px) px[, "Control"])
              Control <- rowMeans(as.data.frame(Control))
              px <- cbind(Control, Case)

              # Calculate 'decision.values' from 'probabilities' using inverse Platt scaling
              px <- as.data.frame(px)
              px <- px[, c("Control", "Case")]
              dv <- as.matrix(log(1 / (1 - px[, "Case"]) - 1))
              colnames(dv) <- "Case/Control"

              # Assign binary values based on probability
              pred <- vector("character", nrow(px))
              pred[px$Case > .5] <- "Case"
              pred[px$Case < .5] <- "Control"

              # Break ties randomly
              case <- sum(array@annot$defineCase == "Case") / nrow(array@annot)
              pred[px$Case == .5] <- sample(c("Case", "Control"),
                                            size = sum(px$Case == .5),
                                            replace = TRUE,
                                            prob = c(case, 1 - case))

              # Clean up pred
              pred <- factor(as.vector(pred), levels = c("Control", "Case"))

              final <- new("ExprsPredict",
                           pred = pred, decision.values = dv, probability = px,
                           actual = array@annot$defineCase)

            }else if(class(array) == "RegrsArray"){

              pred <- colMeans(do.call("rbind", lapply(results, function(x) x@pred)))
              final <- new("RegrsPredict", pred = pred, actual = array@annot$defineCase)

            }else{

              stop("Requested ensemble deployment method not available.")
            }

            if(verbose){

              cat("Ensemble classifier performance:\n")
              print(calcStats(final, aucSkip = TRUE, plotSkip = TRUE))
            }

            return(final)
          }
)
