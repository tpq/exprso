###########################################################
### Duplicate feature selection history

#' Duplicate Feature Selection History
#'
#' Duplicates the feature selection history of a reference \code{ExprsArray}
#'  object. Called directly by \code{predict.ExprsMachine} and indirectly by
#'  \code{predict.ExprsModule}. Ensures that the validation set has undergone
#'  the same steps of feature selection and dimension reduction as the
#'  training set.
#'
#' @param object An \code{ExprsArray} object. The validation set that should undergo
#'  feature selection and dimension reduction.
#' @param reference An \code{ExprsArray} or \code{ExprsModel} object. The object with
#'  feature selection history used as a template for the validation set.
#'
#' @return Returns an \code{ExprsArray} object.
#'
#' @seealso
#' \code{\link{exprso-predict}}
#'
#' @export
setGeneric("modHistory",
           function(object, reference) standardGeneric("modHistory")
)

#' @describeIn modHistory Method to duplicate \code{ExprsArray} feature selection history.
#' @export
setMethod("modHistory", "ExprsArray",
          function(object, reference){

            if(!is.null(object@preFilter)){

              # If reference@preFilter has less (or equal) history than object@preFilter
              if(length(reference@preFilter) <= length(object@preFilter)){

                stop("The provided object does not have less history than the reference object.")
              }

              # If history of reference@preFilter is not also in the object@preFilter
              if(!identical(object@preFilter, reference@preFilter[1:length(object@preFilter)])){

                stop("The provided object history does not match the reference history.")
              }
            }

            # Index first non-overlapping history
            index <- length(object@reductionModel) + 1

            # Manipulate object according to history stored in reference
            for(i in index:length(reference@reductionModel)){

              # If the i-th fs did NOT involve dimension reduction
              if(any(is.na(reference@reductionModel[[i]]))){

                # Build new object
                object <- new(class(object),
                              exprs = object@exprs[reference@preFilter[[i]], , drop = FALSE],
                              annot = object@annot,
                              preFilter = append(object@preFilter,
                                                 list(reference@preFilter[[i]])),
                              reductionModel = append(object@reductionModel,
                                                      list(reference@reductionModel[[i]])))

              }else{

                # Build data according to the i-th feature set
                data <- data.frame(t(object@exprs[reference@preFilter[[i]], ]))

                # Then, apply the i-th reduction model
                if("prcomp" %in% class(reference@reductionModel[[i]])){

                  comps <- predict(reference@reductionModel[[i]], data)

                }else{

                  stop("Uh oh! Reduction model class not recognized.")
                }

                # Build new object
                object <- new(class(object),
                              exprs = t(comps),
                              annot = object@annot,
                              preFilter = append(object@preFilter,
                                                 list(reference@preFilter[[i]])),
                              reductionModel = append(object@reductionModel,
                                                      list(reference@reductionModel[[i]])))
              }
            }

            return(object)
          }
)

###########################################################
### Predict class labels using ExprsModel

#' @name exprso-predict
#' @rdname exprso-predict
#'
#' @title Predict Class Labels
#'
#' @description Predict class labels of a validation set.
#'
#' @details
#'
#' For an \code{ExprsMachine}:
#'
#' An \code{ExprsMachine} object can only predict against an \code{ExprsBinary}
#'  object. An \code{ExprsModule} object can only predict against an
#'  \code{ExprsMulti} object. The validation set should never get modified
#'  once separated from the training set. If the training set used to build
#'  an \code{ExprsModule} had a class missing (i.e., has an NA placeholder),
#'  the \code{ExprsModule} cannot predict the missing class. To learn how
#'  this scenario gets handled, read more at \code{\link{doMulti}}.
#'
#' \code{ExprsPredict} objects store predictions in three slots: \code{@@pred},
#'  \code{@@decision.values}, and \code{@@probability}. The first slot
#'  stores a "final decision" based on the class label with the maximum
#'  predicted probability. The second slot stores a transformation
#'  of the predicted probability for each class calculated by the inverse
#'  of Platt scaling. The predicted probability gets returned by the
#'  \code{predict} method called using the stored \code{@@mach} object.
#'  To learn how these slots get used to calculate classifier performance,
#'  read more at \code{\link{calcStats}}.
#'
#' For an \code{ExprsEnsemble}:
#'
#' At the moment, \code{ExprsEnsemble} can only make predictions
#'  using \code{ExprsMachine} objects. Therefore, it can only predict
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
#'  object (i.e., \code{@@pred}) cast a single (all-or-nothing) vote.
#'  Each subject gets assigned the class that received the most number
#'  of votes (i.e., winner takes all). In both scenarios, ties get
#'  broken randomly with equal weights given to each class.
#'
#' @param object An \code{ExprsModel} object or an \code{ExprsEnsemble} object.
#'  The classifier to deploy.
#' @param array An \code{ExprsArray} object. The validation set.
#' @param how A character string. Select from "probability" or "majority".
#'  See Details for the implication of this choice. Argument applies to
#'  \code{ExprsEnsemble} prediction only.
#' @param verbose A logical scalar. Toggles whether to print from
#'  \code{calcStats}.
#'
#' @return Returns an \code{ExprsPredict} object.
#'
#' @seealso
#' \code{\link{modHistory}}, \code{\link{calcStats}}
NULL

#' @rdname exprso-predict
#' @export
setMethod("predict", "ExprsMachine",
          function(object, array, verbose = TRUE){

            if(class(array) != "ExprsBinary"){

              stop("Uh oh! You can only use an ExprsMachine to predict on an ExprsBinary object.")
            }

            # Reproduce the history stored in the ExprsMachine object
            array <- modHistory(array, reference = object)

            # Build data from modHistory output
            data <- data.frame(t(array@exprs))

            if("naiveBayes" %in% class(object@mach)){

              px <- predict(object@mach, data, type = "raw")

            }else if("lda" %in% class(object@mach)){

              pred <- predict(object@mach, data)
              px <- pred$posterior

            }else if("svm" %in% class(object@mach)){

              pred <- predict(object@mach, data, probability = TRUE)
              px <- attr(pred, "probabilities")

            }else if("nnet" %in% class(object@mach)){

              px <- predict(object@mach, data, type = "raw")
              px <- cbind(px, 1 - px)
              colnames(px) <- c("Case", "Control") # do not delete this line!

            }else if("randomForest" %in% class(object@mach)){

              px <- unclass(predict(object@mach, data, type = "prob"))

            }else if("H2OBinomialModel" %in% class(object@mach)){

              # Import data as H2OFrame via a temporary csv
              colnames(data) <- paste0("id", 1:ncol(data))
              tempFile <- tempfile(fileext = ".csv")
              write.csv(data, tempFile)
              h2o.data <- h2o::h2o.importFile(path = tempFile, destination_frame = "h2o.data")

              # Extract probability
              px <- as.data.frame(h2o::h2o.predict(object@mach, h2o.data))$Case
              px <- cbind(px, 1 - px)
              colnames(px) <- c("Case", "Control") # do not delete this line!

            }else{

              stop("Uh oh! Classification model class not recognized.")
            }

            # Calculate 'decision.values' from 'probability' using inverse Platt scaling
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

            final <- new("ExprsPredict",
                         pred = pred, decision.values = dv, probability = as.matrix(px),
                         actual = array@annot$defineCase)

            if(verbose){

              cat("Individual classifier performance:\n")
              print(calcStats(final, aucSkip = TRUE, plotSkip = TRUE))
            }

            return(final)
          }
)

#' @rdname exprso-predict
#' @export
setMethod("predict", "ExprsModule",
          function(object, array, verbose = TRUE){

            if(class(array) != "ExprsMulti"){

              stop("Uh oh! You can only use an ExprsModule to predict on an ExprsMulti object.")
            }

            if(length(object@mach) != length(levels(array@annot$defineCase))){

              stop("Uh oh! ExprsModule and ExprsMulti must have same number of classes.")
            }

            # For each ExprsMachine stored in @mach
            preds <- vector("list", length(object@mach))
            for(i in 1:length(object@mach)){

              # If the i-th ExprsMachine is missing due to lack of cases during build phase
              if(class(object@mach[[i]]) == "ExprsMachine"){

                # Turn the ExprsMulti object into the i-th ExprsBinary object
                array.i <- array
                array.i@annot$defineCase <- ifelse(as.numeric(array.i@annot$defineCase) == i, "Case", "Control")
                class(array.i) <- "ExprsBinary"

                # Predict the i-th ExprsBinary with the i-th ExprsMachine
                cat("Performing a one-vs-all ExprsMachine prediction with class", i, "set as \"Case\".\n")
                preds[[i]] <- predict(object@mach[[i]], array.i)

              }else if(is.na(object@mach[[i]])){

                cat("The ExprsMachine corresponding to class", i, "is missing. Setting probability to 0.\n")
                missingCase <- rep(0, nrow(array@annot))
                preds[[i]] <- new("ExprsPredict",
                                  pred = as.factor(missingCase),
                                  decision.values = as.matrix(missingCase),
                                  probability = as.matrix(data.frame("Control" = missingCase,
                                                                     "Case" = missingCase)),
                                  actual = array@annot$defineCase)

              }else{

                stop("Uh oh! DEBUG ERROR: 003")
              }
            }

            # Calculate 'decision.values' from 'probability' using inverse Platt scaling
            px <- lapply(preds, function(p) p@probability[, "Case", drop = FALSE])
            dv <- lapply(px, function(prob) log(1 / (1 - prob) - 1))
            px <- do.call(cbind, px); colnames(px) <- 1:ncol(px)
            dv <- do.call(cbind, dv); colnames(dv) <- 1:ncol(dv)

            # calculate weight vector for random sampling during ties
            tieBreaker <- sapply(1:ncol(px), function(case) sum(array@annot$defineCase == case))
            tieBreaker <- tieBreaker / nrow(array@annot)

            # Assign classes based on maximum probability
            pred <- vector("character", nrow(px))
            for(i in 1:nrow(px)){

              # Index maximum probabilities for the i-th subject
              max.i <- which(px[i, ] == max(px[i, ]))

              # If no tie
              if(length(max.i) == 1){

                # Select class with maximum probability
                pred[i] <- as.character(which.max(px[i, ]))

              }else{

                # If none of the tied classes appear in test set, sample of all classes
                if(all(tieBreaker[max.i] %in% 0)) max.i <- 1:ncol(px)

                # Take weighted sample based on weight vector
                pred[i] <- sample(levels(array@annot$defineCase)[max.i], size = 1, prob = tieBreaker[max.i])
              }
            }

            # Clean up pred
            pred <- factor(as.numeric(pred), levels = 1:ncol(px))

            final <- new("ExprsPredict",
                         pred = pred, decision.values = dv, probability = px,
                         actual = array@annot$defineCase)

            if(verbose){

              cat("Multi-class classifier performance:\n")
              print(calcStats(final, aucSkip = TRUE, plotSkip = TRUE))
            }

            return(final)
          }
)

###########################################################
### Calculate classifier performance

#' Calculate Classifier Performance
#'
#' \code{calcStats} calculates classifier performance based on the class predictions
#'  and actual class labels stored in an \code{ExprsPredict} object.
#'
#' This function calculates classifier performance based on the predicted
#'  class labels and the actual class labels in one of two ways. If the argument
#'  \code{aucSkip = FALSE} AND the \code{ExprsArray} object was an \code{ExprsBinary}
#'  object with at least one case and one control AND \code{ExprsPredict} contains
#'  a coherent \code{@@probability} slot, \code{calcStats} will calculate classifier
#'  performance using the area under the receiver operating characteristic (ROC) curve
#'  via the \code{ROCR} package. Otherwise, \code{calcStats} will calculate classifier
#'  performance traditionally using a confusion matrix. Note that accuracies calculated
#'  using \code{ROCR} may differ from accuracies calculated using a confusion
#'  matrix because the former may adjust the discrimination threshold to optimize
#'  sensitivity and specificity. The discrimination threshold is automatically chosen
#'  as the point along the ROC which minimizes the Euclidean distance from (0, 1).
#'
#' @param object An \code{ExprsPredict} object.
#' @param aucSkip A logical scalar. Toggles whether to calculate area under the
#'  receiver operating characteristic curve. See details.
#' @param plotSkip A logical scalar. Toggles whether to plot the receiver
#'  operating characteristic curve.
#'
#' @return Returns a \code{data.frame} of performance metrics.
#'
#' @seealso
#' \code{\link{exprso-predict}}
#'
#' @export
setGeneric("calcStats",
           function(object, aucSkip = FALSE, plotSkip = FALSE) standardGeneric("calcStats")
)

#' @describeIn calcStats Method to calculate the performance of a deployed classifier.
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
