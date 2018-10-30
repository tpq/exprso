#' @name exprso-predict
#' @rdname exprso-predict
#'
#' @title Deploy Model
#'
#' @description Deploy a model to predict outcomes from the data.
#'
#' @details
#' Models can only get deployed on an object of the type used to build
#'  the model. Binary classification and regression are handled natively
#'  by the machine learning algorithm chosen. Multi-class classification
#'  is handled by \code{\link{doMulti}}. Note that a validation set
#'  should never get modified once separated from the training set.
#'  See \code{\link{buildEnsemble}} to learn about ensembles.
#'
#' For binary classifier ensembles, when \code{how = "probability"}, outcomes
#'  are based on the average class probability (via \code{@@probability})
#'  estimated by each deployed model. When \code{how = "majority"}, outcomes
#'  are based on consensus voting whereby each deployed model casts a single
#'  (all-or-nothing) vote (via \code{@@pred}) in a winner takes all approach.
#'  In both scenarios, ties get broken randomly (as weighted by class).
#'
#' For multi-class classifier ensembles, outcomes are based on the
#'  \code{how = "majority"} method from above. For regression ensembles,
#'  outcomes are based on the average predicted value.
#'
#' @param object An \code{ExprsModel} or \code{ExprsEnsemble} object.
#' @param array An \code{ExprsArray} object. The target data.
#' @param how A character string. Select from "probability" or "majority".
#'  See Details. Argument applies to binary classifier ensembles only.
#' @param verbose A logical scalar. Argument passed to \code{calcStats}.
#' @return Returns an \code{ExprsPredict} or \code{RegrsPredict} object.
NULL

#' @rdname exprso-predict
#' @export
setMethod("predict", "ExprsMachine",
          function(object, array, verbose = TRUE){

            classCheck(array, "ExprsBinary",
                       "You must deploy this model on an object built for binary classification.")

            # Reproduce the history stored in the ExprsMachine object
            array <- modHistory(array, reference = object)
            data <- data.frame(t(array@exprs))

            if("naiveBayes" %in% class(object@mach)){

              px <- predict(object@mach, data, type = "raw")

            }else if("lda" %in% class(object@mach)){

              pred <- predict(object@mach, data)
              px <- pred$posterior

            }else if("svm" %in% class(object@mach)){

              pred <- predict(object@mach, data, probability = TRUE)
              px <- attr(pred, "probabilities")

            }else if("lm" %in% class(object@mach)){ # logistic regression

              px <- stats::plogis(predict(object@mach, data))
              px <- cbind(px, 1 - px)
              colnames(px) <- c("Case", "Control") # do not delete this line!

            }else if("cv.glmnet" %in% class(object@mach)){ # LASSO

              px <- stats::plogis(predict(object@mach, as.matrix(data)))
              px <- cbind(px, 1 - px)
              colnames(px) <- c("Case", "Control") # do not delete this line!

            }else if("nnet" %in% class(object@mach)){

              px <- predict(object@mach, data, type = "raw")
              px <- cbind(px, 1 - px)
              colnames(px) <- c("Case", "Control") # do not delete this line!

            }else if("rpart" %in% class(object@mach)){

              px <- predict(object@mach, data)

            }else if("frbs" %in% class(object@mach)){

              px <- round(predict(object@mach, data))
              px[, 1] <- ifelse(px[, 1] == 2, 1, 0) # 1 is Control; 2 is Case
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

              stop("Uh oh! DEBUG ERROR: 004")
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

            classCheck(array, "ExprsMulti",
                       "You must deploy this model on an object built for multi-class classification.")

            if(length(object@mach) != length(levels(array@annot$defineCase))){
              stop("Model and target data not built with same number of classes.")
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

                stop("Uh oh! DEBUG ERROR: 005")
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

#' @rdname exprso-predict
#' @export
setMethod("predict", "RegrsModel",
          function(object, array, verbose = TRUE){

            classCheck(array, "RegrsArray",
                       "You must deploy this model on an object built for continuous outcome tasks.")

            # Reproduce the history stored in the ExprsMachine object
            array <- modHistory(array, reference = object)
            data <- data.frame(t(array@exprs))

            if("svm" %in% class(object@mach) |
               "nnet" %in% class(object@mach) |
               "randomForest" %in% class(object@mach) |
               "rpart" %in% class(object@mach) |
               "lm" %in% class(object@mach)){

              pred <- predict(object@mach, data)

            }else if("cv.glmnet" %in% class(object@mach)){

              pred <- predict(object@mach, as.matrix(data))

            }else if("frbs" %in% class(object@mach)){

              pred <- as.numeric(predict(object@mach, data))

            }else if("H2ORegressionModel" %in% class(object@mach)){

              # Import data as H2OFrame via a temporary csv
              colnames(data) <- paste0("id", 1:ncol(data))
              tempFile <- tempfile(fileext = ".csv")
              write.csv(data, tempFile)
              h2o.data <- h2o::h2o.importFile(path = tempFile, destination_frame = "h2o.data")
              pred <- as.data.frame(h2o::h2o.predict(object@mach, h2o.data))

            }else{

              stop("Uh oh! DEBUG ERROR: 006")
            }

            pred <- as.numeric(as.matrix(pred))
            final <- new("RegrsPredict", pred = pred, actual = array@annot$defineCase)
            if(verbose){
              cat("Continuous outcome performance:\n")
              print(calcStats(final))
            }

            return(final)
          }
)


#' @rdname exprso-predict
#' @export
setMethod("predict", "ExprsEnsemble",
          function(object, array, how = "probability", verbose = TRUE){

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
