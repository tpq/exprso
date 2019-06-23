#' Build Ensemble
#'
#' \code{buildEnsemble} builds an ensemble from \code{ExprsModel} or
#'  \code{ExprsPipeline} objects. See Details.
#'
#' This function can combine any number of model objects into an ensemble.
#'  These models do not necessarily have to derive from the same \code{build}
#'  method. In this way, it works like \code{\link{conjoin}}.
#'
#' This function can also build an ensemble from pipeline objects. It does
#'  this by calling \code{\link{pipeFilter}}, then joining the remaining models
#'  into an ensemble. As an adjunct to this method, consider first combining
#'  multiple pipeline objects with \code{\link{conjoin}}.
#'
#' @inheritParams pipeFilter
#' @param ... Additional \code{ExprsModel} objects to use in the ensemble.
#'  Argument applies to the \code{\link{ExprsModel-class}} method only.
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

            object <- pipeFilter(object, colBy = colBy, how = how, gate = gate, top = top)
            new("ExprsEnsemble", machs = unlist(object@machs))
          }
)

#' @rdname exprso-predict
#'
#' @details
#' For regression ensembles, the average outcome is reported. For multi-class
#'  classifier ensembles, the majority vote is reported. For binary classifier
#'  ensembles, the majority vote or probability-weighted vote is reported.
#'  For probability-weighted voting considers the threshold, the average
#'  "Case" probability is reported. All ties broken randomly.
#'
#' @param how A string. Describes how the ensemble decides. By default, it
#'  uses "majority" voting. However, the user can select "probability" voting
#'  for binary classifier ensembles.
#'
#' @export
setMethod("predict", "ExprsEnsemble",
          function(object, array, how = "majority", verbose = TRUE){

            # Deploy each machine in @machs on the provided ExprsArray
            results <- lapply(object@machs, function(mach) predict(mach, array, verbose = verbose))

            if(class(array) == "ExprsBinary"){

              if(how == "majority"){

                # Majority vote on best outcome
                votes <- data.frame(lapply(results, function(result) result@pred))
                pred <- apply(votes, 1, function(x) names(table(x))[nnet::which.is.max(table(x))])

                # Clean up pred
                pred <- factor(pred, levels = c("Control", "Case"))
                px <- NULL
                dv <- NULL

              }else if(how == "probability"){

                # Get average probabilities
                pxs <- lapply(results, function(result) result@probability)
                Case <- rowMeans(data.frame(lapply(pxs, function(px) px[, "Case"])))
                px <- data.frame("Control" = 1 - Case, "Case" = Case)

                # Get 'decision.values' from 'probabilities' using inverse Platt scaling
                dv <- as.matrix(log(1 / (1 - px[, "Case"]) - 1))
                colnames(dv) <- "Case/Control"

                # Majority vote on best outcome based on average probability
                pred <- ifelse(px$Case > .5, "Case", "Control")
                pred[px$Case == .5] <- sample(c("Case", "Control"))

                # Clean up pred
                pred <- factor(pred, levels = c("Control", "Case"))

              }else{

                stop("Provided 'how' not recognized.")
              }

              final <- new("ExprsPredict",
                           pred = pred, decision.values = dv, probability = px,
                           actual = array$defineCase)

            }else if(class(array) == "ExprsMulti"){

              # Majority vote on best outcome
              votes <- data.frame(lapply(results, function(result) result@pred))
              pred <- apply(votes, 1, function(x) names(table(x))[nnet::which.is.max(table(x))])

              # Clean up pred
              pred <- factor(pred, levels = levels(array$defineCase))
              final <- new("MultiPredict", pred = pred, actual = array$defineCase)

            }else if(class(array) == "RegrsArray"){

              # Average the predictions
              pred <- colMeans(do.call("rbind", lapply(results, function(x) x@pred)))
              final <- new("RegrsPredict", pred = pred, actual = array$defineCase)
            }

            if(verbose){
              cat("Ensemble classifier performance:\n")
              print(calcStats(final, aucSkip = TRUE, plotSkip = TRUE))
            }

            return(final)
          }
)
