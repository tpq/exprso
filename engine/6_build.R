###########################################################
### Define generic functions

setGeneric("buildNB", function(object, ...) standardGeneric("buildNB"))
setGeneric("buildLDA", function(object, ...) standardGeneric("buildLDA"))
setGeneric("buildSVM", function(object, ...) standardGeneric("buildSVM"))
setGeneric("buildANN", function(object, ...) standardGeneric("buildANN"))
setGeneric("buildRF", function(object, ...) standardGeneric("buildRF"))

setGeneric("modHistory", function(object, ...) standardGeneric("modHistory"))
setGeneric("calcStats", function(object, ...) standardGeneric("calcStats"))

###########################################################
### Build classifier

# NOTE: User has one more opportunity to subset the probe set before building machine
# NOTE: @preFilter and @reductionModel data will get passed to ExprsMachine object
# NOTE: Once predict is called, validation set will get processed accordingly

setMethod("buildNB", "ExprsBinary",
          function(object, probes, ...){ # args to naiveBayes
            
            require(e1071)
            
            args <- as.list(substitute(list(...)))[-1]
            
            # Convert 'numeric' probe argument to 'character' probe vector
            if(class(probes) == "numeric"){
              
              if(probes == 0) probes <- nrow(object@exprs)
              probes <- rownames(object@exprs[1:probes, ])
              data <- t(object@exprs[probes, ])
            }
            
            # Build data using supplied 'character' probe vector
            if(class(probes) == "character"){
              
              data <- t(object@exprs[probes, ])
            }
            
            # Set labels as factor
            labels <- factor(object@annot[rownames(data), "defineCase"], levels = c("Control", "Case"))
            
            # Perform naiveBayes via ~ method
            df <- data.frame(data, "defineCase" = labels)
            args <- append(list("formula" = defineCase ~ ., "data" = df), args)
            model <- do.call(naiveBayes, args)
            
            # Carry through and append fs history as stored in the ExprsArray object
            # NOTE: length(ExprsMachine@preFilter) > length(ExprsArray@preFilter)
            machine <- new("ExprsMachine",
                           preFilter = append(object@preFilter, list(probes)),
                           reductionModel = append(object@reductionModel, list(NA)),
                           mach = model
            )
            
            return(machine)
          }
)

setMethod("buildLDA", "ExprsBinary",
          function(object, probes, ...){ # args to lda
            
            require(MASS)
            
            args <- as.list(substitute(list(...)))[-1]
            
            # Convert 'numeric' probe argument to 'character' probe vector
            if(class(probes) == "numeric"){
              
              if(probes == 0) probes <- nrow(object@exprs)
              probes <- rownames(object@exprs[1:probes, ])
              data <- t(object@exprs[probes, ])
            }
            
            # Build data using supplied 'character' probe vector
            if(class(probes) == "character"){
              
              data <- t(object@exprs[probes, ])
            }
            
            # Set labels as factor
            labels <- factor(object@annot[rownames(data), "defineCase"], levels = c("Control", "Case"))
            
            # Perform linear discriminant analysis via ~ method
            df <- data.frame(data, "defineCase" = labels)
            args <- append(list("formula" = defineCase ~ ., "data" = df), args)
            model <- do.call(lda, args)
            
            # Carry through and append fs history as stored in the ExprsArray object
            # NOTE: length(ExprsMachine@preFilter) > length(ExprsArray@preFilter)
            machine <- new("ExprsMachine",
                           preFilter = append(object@preFilter, list(probes)),
                           reductionModel = append(object@reductionModel, list(NA)),
                           mach = model
            )
            
            return(machine)
          }
)

setMethod("buildSVM", "ExprsBinary",
          function(object, probes, ...){ # args to svm
            
            require(e1071)
            
            args <- as.list(substitute(list(...)))[-1]
            
            if(!"probability" %in% names(args)){
              
              cat("Setting 'probability' to TRUE (default behavior, override explicitly)...\n")
              args <- append(args, list("probability" = TRUE))
            }
            
            # Convert 'numeric' probe argument to 'character' probe vector
            if(class(probes) == "numeric"){
              
              if(probes == 0) probes <- nrow(object@exprs)
              probes <- rownames(object@exprs[1:probes, ])
              data <- t(object@exprs[probes, ])
            }
            
            # Build data using supplied 'character' probe vector
            if(class(probes) == "character"){
              
              data <- t(object@exprs[probes, ])
            }
            
            # Set labels as factor
            labels <- factor(object@annot[rownames(data), "defineCase"], levels = c("Control", "Case"))
            
            # Perform SVM via ~ method (permits plotting)
            df <- data.frame(data, "defineCase" = labels)
            args <- append(list("formula" = defineCase ~ ., "data" = df), args)
            model <- do.call(svm, args)
            
            # Carry through and append fs history as stored in the ExprsArray object
            # NOTE: length(ExprsMachine@preFilter) > length(ExprsArray@preFilter)
            machine <- new("ExprsMachine",
                           preFilter = append(object@preFilter, list(probes)),
                           reductionModel = append(object@reductionModel, list(NA)),
                           mach = model
            )
            
            return(machine)
          }
)

# NOTE: Consider varying decay = c(.5, 1) and size = c(1:7)
setMethod("buildANN", "ExprsBinary",
          function(object, probes, ...){
            
            require(nnet)
            
            args <- as.list(substitute(list(...)))[-1]
            
            if(!"size" %in% names(args)){
              
              cat("Setting 'size' to 1 (default behavior, override explicitly)...\n")
              args <- append(args, list("size" = 1))
            }
            
            if(!"range" %in% names(args)){
              
              cat("Setting 'range' to 1/max(|x|) (default behavior, override explicitly)...\n")
              args <- append(args, list("range" = 1/max(abs(as.vector(object@exprs)))))
            }
            
            if(!"decay" %in% names(args)){
              
              cat("Setting 'decay' to 0.5 (default behavior, override explicitly)...\n")
              args <- append(args, list("decay" = 0.5))
            }
            
            if(!"maxit" %in% names(args)){
              
              cat("Setting 'maxit' to 1000 (default behavior, override explicitly)...\n")
              args <- append(args, list("maxit" = 1000))
            }
            
            # Convert 'numeric' probe argument to 'character' probe vector
            if(class(probes) == "numeric"){
              
              if(probes == 0) probes <- nrow(object@exprs)
              probes <- rownames(object@exprs[1:probes, ])
              data <- t(object@exprs[probes, ])
            }
            
            # Build data using supplied 'character' probe vector
            if(class(probes) == "character"){
              
              data <- t(object@exprs[probes, ])
            }
            
            # Set labels as factor
            labels <- factor(object@annot[rownames(data), "defineCase"], levels = c("Control", "Case"))
            
            # Perform ANN via ~ method
            df <- data.frame(data, "defineCase" = labels)
            args <- append(list("formula" = defineCase ~ ., "data" = df), args)
            model <- do.call(nnet, args)
            
            # Carry through and append fs history as stored in the ExprsArray object
            # NOTE: length(ExprsMachine@preFilter) > length(ExprsArray@preFilter)
            machine <- new("ExprsMachine",
                           preFilter = append(object@preFilter, list(probes)),
                           reductionModel = append(object@reductionModel, list(NA)),
                           mach = model
            )
            
            return(machine)
          }
)

setMethod("buildRF", "ExprsBinary",
          function(object, probes, ...){ # args to randomForest
            
            require(randomForest)
            
            args <- as.list(substitute(list(...)))[-1]
            
            # Convert 'numeric' probe argument to 'character' probe vector
            if(class(probes) == "numeric"){
              
              if(probes == 0) probes <- nrow(object@exprs)
              probes <- rownames(object@exprs[1:probes, ])
              data <- t(object@exprs[probes, ])
            }
            
            # Build data using supplied 'character' probe vector
            if(class(probes) == "character"){
              
              data <- t(object@exprs[probes, ])
            }
            
            # Set labels as factor
            labels <- factor(object@annot[rownames(data), "defineCase"], levels = c("Control", "Case"))
            
            # Perform RF via ~ method
            df <- data.frame(data, "defineCase" = labels)
            args <- append(list("formula" = defineCase ~ ., "data" = df), args)
            model <- do.call(randomForest, args)
            
            # Carry through and append fs history as stored in the ExprsArray object
            # NOTE: length(ExprsMachine@preFilter) > length(ExprsArray@preFilter)
            machine <- new("ExprsMachine",
                           preFilter = append(object@preFilter, list(probes)),
                           reductionModel = append(object@reductionModel, list(NA)),
                           mach = model
            )
            
            return(machine)
          }
)

###########################################################
### Predict

# NOTE: The validation set should not get modified once separated from the training set
# NOTE: This function performs only rudimentary checks for this type of error
# NOTE: predict.nnet class predictions match @decision.values, but differ when using ROCR
# NOTE: predict.svm class predictions do not match @decision.values, but match ROCR
setMethod("predict", "ExprsMachine",
          function(object, array, verbose = TRUE, ...){ # args to predict(ExprsMachine, array)
            
            if(class(array) != "ExprsBinary"){
              
              stop("Uh oh! You can only use an ExprsMachine to predict on an ExprsBinary object.")
            }
            
            # Reproduce the history stored in the ExprsMachine object
            array <- modHistory(array, reference = object)
            
            # Build data from modHistory output
            data <- t(array@exprs)
            
            # NOTE: as.data.frame will not rename columns
            data <- data.frame(data)
            
            if("naiveBayes" %in% class(object@mach)){
              
              require(e1071)
              
              # Predict naiveBayes via ~ method
              pred <- predict(object@mach, data, type = "class")
              
              # Calculate 'decision.values' from 'probabilities' using inverse Platt scaling
              px <- predict(object@mach, data, type = "raw")
              px <- as.data.frame(px)
              px <- px[, c("Control", "Case")] # order columns
              dv <- as.matrix(log(1 / (1 - px[, "Case"]) - 1))
              colnames(dv) <- "Case/Control"
              
              # Clean up pred
              pred <- factor(as.vector(pred), levels = c("Control", "Case"))
              names(pred) <- rownames(data)
            }
            
            if("lda" %in% class(object@mach)){
              
              require(MASS)
              
              # Predict linear discriminant analysis via ~ method
              prediction <- predict(object@mach, data)
              
              # Retrieve binary 'pred' class predictions
              pred <- as.character(prediction$class)
              
              # Calculate 'decision.values' from 'probabilities' using inverse Platt scaling
              px <- prediction$posterior
              px <- as.data.frame(px)
              px <- px[, c("Control", "Case")] # order columns
              dv <- as.matrix(log(1 / (1 - px[, "Case"]) - 1))
              colnames(dv) <- "Case/Control"
              
              # Clean up pred
              pred <- factor(as.vector(pred), levels = c("Control", "Case"))
              names(pred) <- rownames(data)
            }
            
            if("svm" %in% class(object@mach)){
              
              require(e1071)
              
              args <- as.list(substitute(list(...)))[-1]
              
              if(!"decision.values" %in% names(args)){
                
                cat("Setting 'decision.values' to TRUE (default behavior, override explicitly)...\n")
                args <- append(args, list("decision.values" = TRUE))
              }
              
              if(!"probability" %in% names(args)){
                
                cat("Setting 'probability' to TRUE (default behavior, override explicitly)...\n")
                args <- append(args, list("probability" = TRUE))
              }
              
              # Predict SVM via ~ method (permits plotting)
              args <- append(list("object" = object@mach, "newdata" = data), args)
              pred <- do.call(predict, args)
              
              # Extract 'probabilities' attribute
              px <- attr(pred, "probabilities")
              px <- as.data.frame(px)
              px <- px[, c("Control", "Case")]
              
              # Calculate 'decision.values' from 'probabilities' using inverse Platt scaling
              if(!is.null(px)){
                
                dv <- as.matrix(log(1 / (1 - px[, "Case"]) - 1))
                colnames(dv) <- "Case/Control"
                
              }else{
                
                dv <- NULL
              }
              
              # Clean up pred
              pred <- factor(as.vector(pred), levels = c("Control", "Case"))
              names(pred) <- rownames(data)
            }
            
            if("nnet" %in% class(object@mach)){
              
              require(nnet)
              
              # Predict ANN via ~ method
              pred <- predict(object@mach, data, type = "class")
              
              # Calculate 'decision.values' from 'probabilities' using inverse Platt scaling
              px <- predict(object@mach, data, type = "raw")
              if(min(px) > .5){ labMin <- ifelse(pred[which.min(px)] == "Case", "Control", "Case")
              }else{ labMin <- pred[which.min(px)] }
              if(max(px) < .5){ labMax <- ifelse(pred[which.max(px)] == "Case", "Control", "Case")
              }else{ labMax <- pred[which.max(px)] }
              if(px[which.min(px)] == .5 & px[which.max(px)] == .5){ labMin <- "Control"; labMax <- "Case" }
              px <- cbind(1 - px, px)
              px <- as.data.frame(px)
              colnames(px) <- c(labMin, labMax)
              px <- px[, c("Control", "Case")] # order columns
              dv <- as.matrix(log(1 / (1 - px[, "Case"]) - 1))
              colnames(dv) <- "Case/Control"
              
              # Clean up pred
              pred <- factor(as.vector(pred), levels = c("Control", "Case"))
              names(pred) <- rownames(data)
            }
            
            if("randomForest" %in% class(object@mach)){
              
              require(randomForest)
              
              # Predict RF via ~ method
              pred <- predict(object@mach, data, type = "response")
              
              # Calculate 'decision.values' from 'probabilities' using inverse Platt scaling
              px <- unclass(predict(object@mach, data, type = "prob"))
              px <- as.data.frame(px)
              px <- px[, c("Control", "Case")] # order columns
              dv <- as.matrix(log(1 / (1 - px[, "Case"]) - 1))
              colnames(dv) <- "Case/Control"
              
              # Clean up pred
              pred <- factor(as.vector(pred), levels = c("Control", "Case"))
              names(pred) <- rownames(data)
            }
            
            final <- new("ExprsPredict", pred = pred, decision.values = dv, probability = as.matrix(px))
            
            if(verbose){
              
              cat("Individual classifier performance:\n")
              cat("(NOTE: This is a preview! Actual performances may differ!)\n")
              print(calcStats(final, array, aucSkip = FALSE, plotSkip = TRUE))
            }
            
            return(final)
          }
)

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
                object <- new("ExprsArray",
                              exprs = as.matrix(object@exprs[reference@preFilter[[i]], ]), # Update @exprs
                              annot = object@annot, # Preserve @annot
                              preFilter = append(object@preFilter, list(reference@preFilter[[i]])), # Append history
                              reductionModel = append(object@reductionModel, list(reference@reductionModel[[i]])))
                
                # Clean up 1-subject artifact
                colnames(object@exprs) <- rownames(object@annot)
                
              }else{
                
                # Build data according to the i-th probe set
                data <- t(object@exprs[reference@preFilter[[i]], ])
                
                # NOTE: as.data.frame will not rename columns
                data <- data.frame(data)
                
                # Then, apply the i-th reduction model
                comps <- predict(reference@reductionModel[[i]], data)
                
                # Build new object
                object <- new(class(object),
                              exprs = t(comps), # Update @exprs
                              annot = object@annot, # Preserve @annot
                              preFilter = append(object@preFilter, list(reference@preFilter[[i]])), # Append history
                              reductionModel = append(object@reductionModel, list(reference@reductionModel[[i]])))
              }
            }
            
            return(object)
          }
)

setMethod("calcStats", "ExprsPredict",
          function(object, array, aucSkip = FALSE, plotSkip = FALSE){
            
            if(!inherits(array, "ExprsArray")){
              
              stop("Uh oh! You can only assess the performance of a classifier using an ExprsArray object.")
            }
            
            layout(matrix(c(1), 1, 1, byrow = TRUE))
            
            # Build 'actual' object as binary
            actual <- as.numeric(array@annot$defineCase == "Case")
            
            # Build 'predicted' object as binary
            predicted <- as.numeric(object@pred == "Case")
            
            # If predicted set contains more than one class and aucSkip = FALSE
            if(length(unique(actual)) > 1 & !aucSkip){
              
              require(ROCR)
              
              # If pred object contains @probability
              if(!is.null(object@probability)){
                
                cat("Calculating accuracy using ROCR based on prediction probabilities...\n")
                p <- prediction(object@probability[, "Case"], actual)
                
                # Plot AUC curve
                perf <- performance(p, measure = "tpr", x.measure = "fpr")
                if(!plotSkip) plot(perf, col = rainbow(10))
                
                # Index optimal cutoff based on Euclidean distance
                index <- which.min(sqrt((1 - perf@y.values[[1]])^2 + (0 - perf@x.values[[1]])^2))
                
                # Calculate performance metrics
                acc <- performance(p, "acc")@y.values[[1]][index]
                sens <- performance(p, "sens")@y.values[[1]][index]
                spec <- performance(p, "spec")@y.values[[1]][index]
                auc <- performance(p, "auc")@y.values[[1]]
                
                return(data.frame(acc, sens, spec, auc))
                
              }else{
                
                cat("Calculating accuracy using ROCR based on raw votes...\n")
                p <- prediction(predicted, actual)
                
                # Plot AUC curve
                perf <- performance(p, measure = "tpr", x.measure = "fpr")
                if(!plotSkip) plot(perf, col = rainbow(10))
                
                # Index optimal cutoff based on Euclidean distance
                index <- which.min(sqrt((1 - perf@y.values[[1]])^2 + (0 - perf@x.values[[1]])^2))
                
                # Calculate performance metrics
                acc <- performance(p, "acc")@y.values[[1]][index]
                sens <- performance(p, "sens")@y.values[[1]][index]
                spec <- performance(p, "spec")@y.values[[1]][index]
                
                return(data.frame(acc, sens, spec))
              }
              
            }else{
              
              if(!aucSkip){ cat("Predictions include only one class. Calculating accuracy outside of ROCR...\n")
              }else{ cat("aucSkip set to TRUE. Calculating accuracy outside of ROCR...\n")}
              table <- matrix(0, nrow = 2, ncol = 2)
              
              # Fill in the 2 x 2 table
              for(i in 1:nrow(array@annot)){
                
                if(predicted[i] == 1 & actual[i] == 1) table[1, 1] <- table[1, 1] + 1
                if(predicted[i] == 1 & actual[i] == 0) table[1, 2] <- table[1, 2] + 1
                if(predicted[i] == 0 & actual[i] == 1) table[2, 1] <- table[2, 1] + 1
                if(predicted[i] == 0 & actual[i] == 0) table[2, 2] <- table[2, 2] + 1
              }
              
              # Calculate performance metrics
              acc <- (table[1, 1] + table[2, 2]) / sum(table)
              sens <- (table[1, 1]) / (table[1, 1] + table[2, 1])
              spec <- (table[2, 2]) / (table[1, 2] + table[2, 2])
              
              return(data.frame(acc, sens, spec))
            }
          }
)
