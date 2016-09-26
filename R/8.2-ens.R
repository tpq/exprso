###########################################################
### Build Ensemble

#' Build Ensemble
#'
#' Aggregates multiple classifiers into a single ensemble classifier.
#'
#' The \code{\link{ExprsModel-class}} method:
#'
#' Combine any number of \code{ExprsModel} objects into an ensemble. These models
#'  do not necessarily have to derive from the same \code{build} method. This
#'  method works identically to the \code{\link{conjoin}} \code{ExprsModel} method.
#'
#' The \code{\link{ExprsPipeline-class}} method:
#'
#' Build an ensemble from an \code{ExprsPipeline} object. This method works by
#'  calling \code{\link{pipeFilter}}, then aggregating those results into an ensemble.
#'  As an adjunct to this method, consider first combining multiple
#'  \code{ExprsPipeline} objects together with \code{\link{conjoin}}.
#'
#' @param object An \code{\link{ExprsModel-class}} object.
#' @param ... Additional \code{ExprsModel} objects to use in the ensemble.
#'  Argument applies to \code{\link{ExprsModel-class}} method only.
#' @inheritParams pipeFilter
#'
#' @return An \code{\link{ExprsEnsemble-class}} object.
#'
#' @seealso
#' \code{\link{pipeFilter}}\cr
#' \code{\link{pipeUnboot}}\cr
#' \code{\link{plCV}}\cr
#' \code{\link{plGrid}}\cr
#' \code{\link{plGridMulti}}\cr
#' \code{\link{plMonteCarlo}}\cr
#' \code{\link{plNested}}
#'
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

###########################################################
### Predict class labels using ExprsEnsemble

#' @rdname exprso-predict
#' @export
setMethod("predict", "ExprsEnsemble",
          function(object, array, how = "probability", verbose = TRUE){

            if(!inherits(array, "ExprsArray")){

              stop("Uh oh! You can only use an ExprsEnsemble to predict on an ExprsArray object!")
            }

            # Deploy each machine in @machs on the provided ExprsArray
            results <- lapply(object@machs, function(mach) predict(mach, array, verbose = verbose))

            if("ExprsBinary" %in% class(array)){

              if(how == "probability"){

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

              }else if(how == "majority"){

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

              }else{

                stop("Uh oh! Provided 'how' argument not recognized!")
              }

            }else if("ExprsMulti" %in% class(array)){

              stop("Uh oh! predict.ExprsEnsemble not yet generalized to work for ExprsMulti objects!")

            }else{

              stop("Uh oh! You can only use an ExprsEnsemble to predict on an ExprsArray object!")
            }

            final <- new("ExprsPredict",
                         pred = pred, decision.values = dv, probability = px,
                         actual = array@annot$defineCase)

            if(verbose){

              cat("Ensemble classifier performance:\n")
              print(calcStats(final, aucSkip = TRUE, plotSkip = TRUE))
            }

            return(final)
          }
)
