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
#' For regression, accuracy is defined the R-squared of the fitted regression. This
#'  ranges from 0 to 1 for use with \code{\link{pl}} and \code{\link{pipe}}. Note that
#'  the \code{aucSkip} and \code{plotSkip} arguments are ignored for regression.
#'
#' @param object An \code{ExprsPredict} or \code{RegrsPredict} object.
#' @param aucSkip A logical scalar. Toggles whether to calculate area under the
#'  receiver operating characteristic curve. See Details.
#' @param plotSkip A logical scalar. Toggles whether to plot the receiver
#'  operating characteristic curve. See Details.
#' @param verbose A logical scalar. Toggles whether to print the results
#'  of model performance to console.
#'
#' @return Returns a \code{data.frame} of performance metrics.
#'
#' @export
setGeneric("calcStats",
           function(object, aucSkip = FALSE, plotSkip = FALSE, verbose = TRUE) standardGeneric("calcStats")
)

#' @describeIn calcStats Method to calculate performance for classification models.
#' @export
setMethod("calcStats", "ExprsPredict",
          function(object, aucSkip, plotSkip, verbose){

            if(all(c("Case", "Control") %in% object@actual) & !is.null(object@probability) & !aucSkip){

              if(verbose) cat("Calculating accuracy based on optimal AUC cutoff...\n")

              # Find optimal cutoff based on distance from top-left corner
              p <- ROCR::prediction(object@probability[, "Case"], as.numeric(object@actual == "Case"))
              perf <- ROCR::performance(p, measure = "tpr", x.measure = "fpr")
              index <- which.min(sqrt((1 - perf@y.values[[1]])^2 + (0 - perf@x.values[[1]])^2))
              if(!plotSkip) plot(perf, col = rainbow(10))
              if(!plotSkip) graphics::points(perf@x.values[[1]][index], perf@y.values[[1]][index], col = "blue")

              # Get performance for optimal cutoff
              acc <- ROCR::performance(p, "acc")@y.values[[1]][index]
              sens <- ROCR::performance(p, "sens")@y.values[[1]][index]
              spec <- ROCR::performance(p, "spec")@y.values[[1]][index]
              prec <- ROCR::performance(p, "prec")@y.values[[1]][index]
              f1 <- 2 * (prec * sens) / (prec + sens)
              auc <- ROCR::performance(p, "auc")@y.values[[1]]

              df <- data.frame(acc, sens, spec, prec, f1, auc)
              df[is.na(df)] <- 0
              return(df)

            }else{

              if(verbose) cat("Calculating accuracy without AUC support...\n")

              # Build confusion table from factor
              object@actual <- factor(object@actual, levels = c("Control", "Case"))
              table <- table("predicted" = object@pred, "actual" = object@actual)
              if(verbose){
                cat("Classification confusion table:\n"); print(table)
              }

              tn <- table[1,1]
              fp <- table[2,1]
              fn <- table[1,2]
              tp <- table[2,2]

              acc <- (tp + tn) / (tp + tn + fp + fn)
              sens <- tp / (tp + fn)
              spec <- tn / (fp + tn)
              prec <- tp / (tp + fp)
              f1 <- 2 * (prec * sens) / (prec + sens)

              df <- data.frame(acc, sens, spec, prec, f1)
              df[is.na(df)] <- 0

              if(verbose){
                cat("Classification confusion table:\n"); print(table)
                cat("Classifier model performance:\n"); print(df)
              }

              return(df)
            }
          }
)


#' @describeIn calcStats Method to calculate performance for multi-class models.
#' @export
setMethod("calcStats", "MultiPredict",
          function(object, verbose){

            mat <- table("predicted" = object@pred, "actual" = object@actual) # predicted as rows
            acc <- sum(diag(mat)) / sum(mat)
            prec <- diag(mat) / rowSums(mat)
            sens <- diag(mat) / colSums(mat)
            df <- data.frame(acc, "sens" = t(sens), "prec" = t(prec)) # ensures 1 row
            df[is.na(df)] <- 0

            if(verbose){
              cat("Classification confusion table:\n"); print(mat)
              cat("Classifier model performance:\n"); print(df)
            }

            return(df)
          }
)

#' @describeIn calcStats Method to calculate performance for continuous outcome models.
#' @export
setMethod("calcStats", "RegrsPredict",
          function(object, verbose){

            mse <- mean((object@pred - object@actual)^2)
            rmse <- sqrt(mse)
            mae <- mean(abs(object@pred - object@actual))
            cor <- stats::cor(object@pred, object@actual, method = "pearson")
            R2 <- cor^2
            acc <- R2
            df <- data.frame(acc, mse, rmse, mae, cor, R2)
            df[is.na(df)] <- 0

            if(verbose){
              cat("Regression model performance:\n"); print(df)
            }

            return(df)
          }
)
