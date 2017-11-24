#' Calculate Model Performance
#'
#' \code{calcStats} calculates the performance of a deployed model.
#'
#' For classification, if the argument \code{aucSkip = FALSE} AND the \code{ExprsArray}
#'  object was an \code{ExprsBinary} object with at least one case and one control AND
#'  \code{ExprsPredict} contains a coherent \code{@@probability} slot, \code{calcStats}
#'  will calculate classifier performance using the area under the receiver operating
#'  characteristic (ROC) curve via the \code{ROCR} package. Otherwise, \code{calcStats}
#'  will calculate classifier performance traditionally using a confusion matrix.
#'  Note that accuracies calculated using \code{ROCR} may differ from those calculated
#'  using a confusion matrix because \code{ROCR} adjusts the discrimination threshold to
#'  optimize sensitivity and specificity. This threshold is automatically chosen as the
#'  point along the ROC which minimizes the Euclidean distance from (0, 1).
#'
#' For regression, accuracy is defined as (1 - MSE / (1 + MSE)). This allows MSE to
#'  range from 0 to 1 for use with \code{\link{pl}} and \code{\link{pipe}}. Note that
#'  the \code{aucSkip} and \code{plotSkip} arguments are ignored for regression.
#'
#' @param object An \code{ExprsPredict} or \code{RegrsPredict} object.
#' @param aucSkip A logical scalar. Toggles whether to calculate area under the
#'  receiver operating characteristic curve. See Details.
#' @param plotSkip A logical scalar. Toggles whether to plot the receiver
#'  operating characteristic curve. See Details.
#' @return Returns a \code{data.frame} of performance metrics.
#' @export
setGeneric("calcStats",
           function(object, aucSkip = FALSE, plotSkip = FALSE) standardGeneric("calcStats")
)

#' @describeIn calcStats Method to calculate performance for classification models.
#' @export
setMethod("calcStats", "ExprsPredict",
          function(object, aucSkip, plotSkip){

            # If predicted set contains only two classes, ExprsPredict has @probability, and aucSkip = FALSE
            if(all(c("Case", "Control") %in% object@actual) & !is.null(object@probability) & !aucSkip){

              cat("Calculating accuracy using ROCR based on prediction probabilities...\n")

              p <- ROCR::prediction(object@probability[, "Case"], as.numeric(object@actual == "Case"))

              # Plot AUC curve and index optimal cutoff based on Euclidean distance
              perf <- ROCR::performance(p, measure = "tpr", x.measure = "fpr")
              if(!plotSkip) plot(perf, col = rainbow(10))
              index <- which.min(sqrt((1 - perf@y.values[[1]])^2 + (0 - perf@x.values[[1]])^2))

              # Calculate performance metrics
              acc <- ROCR::performance(p, "acc")@y.values[[1]][index]
              sens <- ROCR::performance(p, "sens")@y.values[[1]][index]
              spec <- ROCR::performance(p, "spec")@y.values[[1]][index]
              auc <- ROCR::performance(p, "auc")@y.values[[1]]

              df <- data.frame(acc, sens, spec, auc)
              df[is.na(df)] <- 0
              return(df)

            }else{

              cat("Arguments not provided in an ROCR AUC format. Calculating accuracy outside of ROCR...\n")

              # Turn ExprsBinary $defineCase into factor
              if(any(c("Control", "Case") %in% object@actual)){

                object@actual <- factor(object@actual, levels = c("Control", "Case"))
              }

              # Build confusion table
              table <- table("predicted" = object@pred, "actual" = object@actual)
              cat("Classification confusion table:\n"); print(table)

              # Compute per-class performance
              for(class in 1:nrow(table)){

                tp <- sum(table[row(table) == col(table)][class])
                tn <- sum(table[row(table) == col(table)][-class])
                fp <- sum(table[row(table) == class & col(table) != class]) # called class but not
                fn <- sum(table[row(table) != class & col(table) == class]) # is but not called
                acc <- (tp + tn) / (tp + tn + fp + fn)
                sens <- tp / (tp + fn)
                spec <- tn / (fp + tn)

                if(length(levels(object@actual)) > 2){

                  cat("Class", class, "performance (acc, sens, spec):", paste0(acc,", ",sens,", ", spec), "\n")

                }else{

                  # NOTE: class == 2 refers to "Case"
                  if(class == 2){
                    df <- data.frame(acc, sens, spec)
                    df[is.na(df)] <- 0
                    return(df)
                  }
                }
              }

              # Compute total accuracy
              tp <- sum(table[row(table) == col(table)])
              tn <- sum(table[row(table) == col(table)])
              fp <- sum(table[row(table) != col(table)])
              fn <- sum(table[row(table) != col(table)])
              acc <- (tp + tn) / (tp + tn + fp + fn)

              cat("Total accuracy of ExprsModule:", acc, "\n")

              df <- data.frame(acc)
              df[is.na(df)] <- 0
              return(df)
            }
          }
)

#' @describeIn calcStats Method to calculate performance for continuous outcome models.
#' @export
setMethod("calcStats", "RegrsPredict",
          function(object, aucSkip, plotSkip){

            mse <- mean((object@pred - object@actual)^2)
            rmse <- sqrt(mse)
            mae <- mean(abs(object@pred - object@actual))
            acc <- 1 - mse / (1 + mse) # default for pl calls
            df <- data.frame(acc, mse, rmse, mae)
            return(df)
          }
)
