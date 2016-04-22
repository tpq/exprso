###########################################################
### Build Ensemble

#' Build Ensemble
#'
#' Aggregates multiple classifiers into a single ensemble classifier.
#'
#' The \code{buildEnsemble} function allows a user to quickly build an ensemble
#'  classifier. Using the \code{\link{ExprsModel-class}} method, one can specify
#'  manually any number of \code{ExprsModel} objects to include in the ensemble.
#'  These models do not necessarily have to derive from the same \code{build}
#'  method.
#'
#' Using the \code{\link{ExprsPipeline-class}} method, one can also quickly
#'  build an ensemble from an \code{ExprsPipeline} object. This method works
#'  by first passing all arguments to \code{\link{pipeFilter}}, then aggregating
#'  the remaining models into a classifier. As an adjunct to this method, consider
#'  joining multiple \code{ExprsPipeline} objects with \code{\link{conjoin}}.
#'
#' @seealso
#' \code{\link{plCV}}, \code{\link{plGrid}}, \code{\link{plMonteCarlo}}, \code{\link{plNested}}
#'
#' @export
setGeneric("buildEnsemble",
           function(object, ...) standardGeneric("buildEnsemble")
)

#' @describeIn buildEnsemble Method to build ensemble from \code{ExprsModel} objects.
#' @param object An \code{\link{ExprsModel-class}} object.
#' @param ... Additional \code{ExprsModel} objects to use in an ensemble.
#' @export
setMethod("buildEnsemble", "ExprsModel",
          function(object, ...){ # args to include additional ExprsMachine objects

            conjoin(object, ...)
          }
)

#' @describeIn buildEnsemble Method to build ensemble from \code{ExprsPipeline} objects.
#' @inheritParams pipeFilter
#' @export
setMethod("buildEnsemble", "ExprsPipeline",
          function(object, colBy = 0, how = 0, gate = 0, top.N = 0){

            if(!identical(colBy, 0)){

              object <- pipeFilter(object, colBy = colBy, how = how, gate = gate, top.N = top.N)
            }

            new("ExprsEnsemble",
                machs = unlist(object@machs)
            )
          }
)

###########################################################
### Predict class labels using ExprsEnsemble

#' @rdname exprso-predict
#'
#' @details
#'
#' At the moment, \code{ExprsEnsemble} can only make predictions
#'  with \code{ExprsMachine} objects. Therefore, it can only predict
#'  against \code{ExprsBinary} objets. Predicting with ensembles poses
#'  a unique challenge with regard to how to translate multiple
#'  performance scores (one for each classifier in the ensemble) into
#'  a single performance score (for the ensemble as a whole). For now,
#'  the \code{ExprsEnsemble} \code{predict} method offers two options,
#'  toggled with the argument \code{how}. Regardless of the chosen
#'  \code{how}, \code{buildEnsemble} begins by deploying each constituent
#'  classifier on the validation set to yield a list of \code{ExprsPredict}
#'  objects.
#'
#' When \code{how = "probability"}, this method will take the average
#'  predicted class probability (i.e., \code{@@probability} for each
#'  returned \code{ExprsPredict} object (corresponding to each constituent
#'  \code{ExprsModel} object). When \code{how = "majority"}, this method
#'  will let the final decision from each returned \code{ExprsPredict}
#'  object (i.e., \code{@@binary}) cast a single (all-or-nothing) vote.
#'  Each subject gets assigned the class that received the most number
#'  of votes (i.e., winner takes all). In both scenarios, ties get
#'  broken randomly with equal weights given to each class.
#'
#' @param object An \code{ExprsEnsemble} object. The classifier to deploy.
#' @param how A character string. Select from "probability" or "majority".
#'  See "Details" for the implication of this choice.
#'
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
                names(pred) <- rownames(data)

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

            final <- new("ExprsPredict", pred = pred, decision.values = dv, probability = px)

            if(verbose){

              cat("Ensemble classifier performance:\n")
              print(calcStats(final, array, aucSkip = TRUE, plotSkip = TRUE))
            }

            return(final)
          }
)
