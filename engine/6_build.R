###########################################################
### Define generic functions

setGeneric("buildNB", function(object, ...) standardGeneric("buildNB"))
setGeneric("buildLDA", function(object, ...) standardGeneric("buildLDA"))
setGeneric("buildSVM", function(object, ...) standardGeneric("buildSVM"))
setGeneric("buildANN", function(object, ...) standardGeneric("buildANN"))
setGeneric("buildRF", function(object, ...) standardGeneric("buildRF"))

setGeneric("doMulti", function(object, ...) standardGeneric("doMulti"))
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
              
              args <- append(args, list("probability" = TRUE))
            }
            
            # Mandate probability = TRUE
            if(args$probability == FALSE){
              
              cat("Uh oh! This function requires 'probability' = TRUE. Setting 'probability' to TRUE...\n")
              args$probability <- TRUE
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

# doMulti depends on the initial number of levels in the $defineCase factor.
# If an ExprsMulti subset (e.g. from cross-validation) lacks one (or more) cases
# entirely, that ExprsMachine gets replaced with a NULL place-holder.
setMethod("doMulti", "ExprsMulti",
          function(object, probes, what, ...){
            
            args <- as.list(substitute(list(...)))[-1]
            
            # Initialize multi container
            multi <- vector("list", length(levels(object@annot$defineCase)))
            
            # Perform N binary tasks
            for(i in 1:length(levels(object@annot$defineCase))){
              
              # If the i-th ExprsMachine would not have any representative cases
              if(all(!as.numeric(object@annot$defineCase) == i)){
                
                cat(paste0("Missing a representative of class ", i, ". Object replaced with NULL placeholder.\n"))
                multi[[i]] <- NULL
                
              }else{
                
                # Turn the ExprsMulti object into the i-th ExprsBinary object
                temp <- object
                temp@annot$defineCase <- ifelse(as.numeric(temp@annot$defineCase) == i, "Case", "Control")
                class(temp) <- "ExprsBinary"
                
                # Perform the binary task
                cat("Performing a one-vs-all binary task with class", i, "set as \"Case\".\n")
                args.i <- append(list("object" = temp, "probes" = probes), args)
                multi[[i]] <- do.call(what = what, args = args.i)
              }
            }
            
            return(multi)
          }
)

setMethod("buildNB", "ExprsMulti",
          function(object, probes, ...){
            
            # Pass arguments to doMulti
            machs <- doMulti(object, probes, what = "buildNB", ...)
            
            # Carry through and append fs history as stored in the ExprsArray object
            new("ExprsModule",
                preFilter = append(object@preFilter, list(probes)),
                reductionModel = append(object@reductionModel, list(NA)),
                mach = machs)
          }
)

setMethod("buildLDA", "ExprsMulti",
          function(object, probes, ...){
            
            # Pass arguments to doMulti
            machs <- doMulti(object, probes, what = "buildLDA", ...)
            
            # Carry through and append fs history as stored in the ExprsArray object
            new("ExprsModule",
                preFilter = append(object@preFilter, list(probes)),
                reductionModel = append(object@reductionModel, list(NA)),
                mach = machs)
          }
)

setMethod("buildSVM", "ExprsMulti",
          function(object, probes, ...){
            
            # Pass arguments to doMulti
            machs <- doMulti(object, probes, what = "buildSVM", ...)
            
            # Carry through and append fs history as stored in the ExprsArray object
            new("ExprsModule",
                preFilter = append(object@preFilter, list(probes)),
                reductionModel = append(object@reductionModel, list(NA)),
                mach = machs)
          }
)

setMethod("buildANN", "ExprsMulti",
          function(object, probes, ...){
            
            # Pass arguments to doMulti
            machs <- doMulti(object, probes, what = "buildANN", ...)
            
            # Carry through and append fs history as stored in the ExprsArray object
            new("ExprsModule",
                preFilter = append(object@preFilter, list(probes)),
                reductionModel = append(object@reductionModel, list(NA)),
                mach = machs)
          }
)

setMethod("buildRF", "ExprsMulti",
          function(object, probes, ...){
            
            # Pass arguments to doMulti
            machs <- doMulti(object, probes, what = "buildRF", ...)
            
            # Carry through and append fs history as stored in the ExprsArray object
            new("ExprsModule",
                preFilter = append(object@preFilter, list(probes)),
                reductionModel = append(object@reductionModel, list(NA)),
                mach = machs)
          }
)

###########################################################
### Predict

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
                              exprs = object@exprs[reference@preFilter[[i]], , drop = FALSE], # Update @exprs
                              annot = object@annot, # Preserve @annot
                              preFilter = append(object@preFilter, list(reference@preFilter[[i]])), # Append history
                              reductionModel = append(object@reductionModel, list(reference@reductionModel[[i]])))
                
              }else{
                
                # Build data according to the i-th probe set
                data <- data.frame(t(object@exprs[reference@preFilter[[i]], ]))
                
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

# NOTE: The validation set should not get modified once separated from the training set
# NOTE: All @pred and @decision.values now based on @probability
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
              
              require(e1071)
              px <- predict(object@mach, data, type = "raw")
            }
            
            if("lda" %in% class(object@mach)){
              
              require(MASS)
              pred <- predict(object@mach, data)
              px <- pred$posterior
            }
            
            if("svm" %in% class(object@mach)){
              
              require(e1071)
              pred <- predict(object@mach, data, probability = TRUE)
              px <- attr(pred, "probabilities")
              
            }
            
            if("nnet" %in% class(object@mach)){
              
              require(nnet)
              px <- predict(object@mach, data, type = "raw")
              px <- cbind(px, 1 - px)
              colnames(px) <- c("Case", "Control") # do not delete this line!
            }
            
            if("randomForest" %in% class(object@mach)){
              
              require(randomForest)
              px <- unclass(predict(object@mach, data, type = "prob"))
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
            pred[px$Case == .5] <- sample(c("Case", "Control"), sum(px$Case == .5), TRUE, prob = c(case, 1 - case))
            
            # Clean up pred
            pred <- factor(as.vector(pred), levels = c("Control", "Case"))
            names(pred) <- rownames(data)
            
            final <- new("ExprsPredict", pred = pred, decision.values = dv, probability = as.matrix(px))
            
            if(verbose){
              
              cat("Individual classifier performance:\n")
              print(calcStats(final, array, aucSkip = TRUE, plotSkip = TRUE))
            }
            
            return(final)
          }
)

# If the training set used to build the ExprsModule had a class missing (e.g. from cross-validation),
# the ExprsModule will not predict that class. As with all functions included in this package,
# relative class frequencies determine the 'tieBreaker' probability weights.
setMethod("predict", "ExprsModule",
          function(object, array, verbose = TRUE){
            
            if(class(array) != "ExprsMulti"){
              
              stop("Uh oh! You can only use an ExprsModule to predict on an ExprsMulti object.")
            }
            
            if(length(object@mach) != length(levels(array@annot$defineCase))){
              
              stop("Uh oh! ExprsModule and ExprsMulti must have same number of classes.")
            }
            
            # Initialize preds container
            preds <- vector("list", length(object@mach))
            
            # For each ExprsMachine stored in @mach
            for(i in 1:length(object@mach)){
              
              # If the i-th ExprsMachine is missing due to lack of cases during build phase
              if(is.null(object@mach[[i]])){
                
                cat("The ExprsMachine corresponding to class", i, "is missing. Setting probability to 0.\n")
                missingCase <- rep(0, nrow(array@annot))
                preds[[i]] <- new("ExprsPredict", pred = as.factor(missingCase), decision.values = as.matrix(missingCase),
                                  probability = as.matrix(data.frame("Control" = missingCase, "Case" = missingCase)))
              }else{
                
                # Turn the ExprsMulti object into the i-th ExprsBinary object
                array.i <- array
                array.i@annot$defineCase <- ifelse(as.numeric(array.i@annot$defineCase) == i, "Case", "Control")
                class(array.i) <- "ExprsBinary"
                
                # Predict the i-th ExprsBinary with the i-th ExprsMachine
                cat("Performing a one-vs-all ExprsMachine prediction with class", i, "set as \"Case\".\n")
                preds[[i]] <- predict(object@mach[[i]], array.i)
              }
            }
            
            # Calculate 'decision.values' from 'probability' using inverse Platt scaling
            px <- lapply(preds, function(p) p@probability[, "Case", drop = FALSE])
            dv <- lapply(px, function(prob) log(1 / (1 - prob) - 1))
            px <- do.call(cbind, px); colnames(px) <- 1:ncol(px)
            dv <- do.call(cbind, dv); colnames(dv) <- 1:ncol(dv)
            
            # Initialize pred container
            pred <- vector("character", nrow(px))
            
            # calculate weight vector for random sampling during ties
            tieBreaker <- sapply(1:ncol(px), function(case) sum(array@annot$defineCase == case)) / nrow(array@annot)
            
            # Assign classes based on maximum probability
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
                pred[i] <- sample(levels(array@annot$defineCase)[max.i], prob = tieBreaker[max.i])
              }
            }
            
            # Clean up pred
            pred <- factor(as.numeric(pred), levels = 1:ncol(px))
            
            final <- new("ExprsPredict", pred = pred, decision.values = dv, probability = px)
            
            if(verbose){
              
              cat("Multi-class classifier performance:\n")
              print(calcStats(final, array, aucSkip = TRUE, plotSkip = TRUE))
            }
            
            return(final)
          }
)

setMethod("calcStats", "ExprsPredict",
          function(object, array, aucSkip = FALSE, plotSkip = FALSE){
            
            if(!inherits(array, "ExprsArray")){
              
              stop("Uh oh! You can only assess the performance of a classifier using an ExprsArray object.")
            }
            
            if(class(array) == "ExprsMulti"){
              
              if(length(levels(object@pred)) != length(levels(array@annot$defineCase))){
                
                stop("Uh oh! ExprsPredict and ExprsMulti must have same number of classes.")
              }
            }
            
            # If predicted set contains only two classes, ExprsPredict has @probability, and aucSkip = FALSE
            if(all(c("Case", "Control") %in% array@annot$defineCase) & !is.null(object@probability) & !aucSkip){
              
              require(ROCR)
              
              layout(matrix(c(1), 1, 1, byrow = TRUE))
              
              cat("Calculating accuracy using ROCR based on prediction probabilities...\n")
              p <- prediction(object@probability[, "Case"], as.numeric(array@annot$defineCase == "Case"))
              
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
              
              cat("Arguments not provided in an ROCR AUC format. Calculating accuracy outside of ROCR...\n")
              
              # Turn ExprsBinary $defineCase into factor
              if(class(array) == "ExprsBinary"){
                
                array@annot$defineCase <- factor(array@annot$defineCase, levels = c("Control", "Case"))
              }
              
              # Build confusion table
              table <- table("predicted" = object@pred, "actual" = array@annot$defineCase)
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
                
                # If multi-class
                if(class(array) == "ExprsMulti"){
                  
                  cat("Class", class, "performance (acc, sens, spec):", paste0(acc,", ",sens,", ", spec), "\n")
                  
                }else{
                  
                  # NOTE: class == 2 refers to "Case"
                  if(class == 2) return(data.frame(acc, sens, spec))
                }
              }
              
              # Compute total accuracy
              tp <- sum(table[row(table) == col(table)])
              tn <- sum(table[row(table) == col(table)])
              fp <- sum(table[row(table) != col(table)])
              fn <- sum(table[row(table) != col(table)])
              acc <- (tp + tn) / (tp + tn + fp + fn)
              
              cat("Total accuracy of ExprsModule:", acc, "\n")
              
              return(data.frame(acc))
            }
          }
)
