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
              final <- new("MultiPredict", pred = pred, actual = array$defineCase)

              # Clean up pred
              levels <- levels(array$defineCase)
              pred <- factor(pred, levels = 1:length(levels), labels = levels)

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
