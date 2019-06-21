#' @name exprso-predict
#' @rdname exprso-predict
#'
#' @title Deploy Model
#'
#' @description Deploy a model to predict outcomes from the data.
#'
#' @details
#' Models can only get deployed on an object of the type used to build
#'  the model. This method now supports binary classification,
#'  multi-class classification, and regression.
#'
#' @param object An \code{exprso} model.
#' @param array An \code{exprso} object. The test data.
#' @param verbose A boolean. Argument passed to \code{calcStats}.
#' @return Returns an \code{exprso} prediction object.
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

              stop("Uh oh! DEBUG ERROR: PRED1")
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

            # Reproduce the history stored in the ExprsMachine object
            array <- modHistory(array, reference = object)
            data <- data.frame(t(array@exprs))

            if("naiveBayes" %in% class(object@mach) |
               "svm" %in% class(object@mach) |
               "randomForest" %in% class(object@mach)){

              pred <- predict(object@mach, data)

            }else if("lda" %in% class(object@mach)){

              pred <- predict(object@mach, data)
              px <- pred$posterior
              pred <- apply(px, 1, which.max)
              levels <- levels(array$defineCase)
              pred <- factor(pred, levels = 1:length(levels), labels = levels)

            }else if("cv.glmnet" %in% class(object@mach) |
                     "nnet" %in% class(object@mach)){

              # These all return a probability matrix; convert to factor label
              px <- predict(object@mach, as.matrix(data)) # needs matrix
              pred <- apply(px, 1, which.max)
              levels <- levels(array$defineCase)
              pred <- factor(pred, levels = 1:length(levels), labels = levels)

            }else if("rpart" %in% class(object@mach)){

              # These all return a probability matrix; convert to factor label
              px <- predict(object@mach, data)
              pred <- apply(px, 1, which.max)
              levels <- levels(array$defineCase)
              pred <- factor(pred, levels = 1:length(levels), labels = levels)

            }else if("frbs" %in% class(object@mach)){

              # FRBS returns a number; convert to factor label
              pred <- round(predict(object@mach, data))
              levels <- levels(array$defineCase)
              pred <- factor(pred, levels = 1:length(levels), labels = levels)

            }else if("H2ORegressionModel" %in% class(object@mach)){

              # Import data as H2OFrame via a temporary csv
              colnames(data) <- paste0("id", 1:ncol(data))
              tempFile <- tempfile(fileext = ".csv")
              write.csv(data, tempFile)
              h2o.data <- h2o::h2o.importFile(path = tempFile, destination_frame = "h2o.data")
              pred <- as.data.frame(h2o::h2o.predict(object@mach, h2o.data))$predict

            }else{

              stop("Uh oh! DEBUG ERROR: PRED2")
            }

            names(pred) <- NULL
            final <- new("MultiPredict", pred = pred, actual = array@annot$defineCase)
            if(verbose){
              cat("Multi-class classifier performance:\n")
              print(calcStats(final))
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

              stop("Uh oh! DEBUG ERROR: PRED3")
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
