###########################################################
### Define generic functions

setGeneric("pipeUnboot", function(object, ...) standardGeneric("pipeUnboot"))
setGeneric("pipeSubset", function(object, ...) standardGeneric("pipeSubset"))
setGeneric("pipeFilter", function(object, ...) standardGeneric("pipeFilter"))

setGeneric("buildEnsemble", function(object, ...) standardGeneric("buildEnsemble"))

###########################################################
### Modify ExprsPipeline objects

setMethod("pipeUnboot", "ExprsPipeline",
          function(object){
            
            if("boot" %in% colnames(object@summary)){
              
              # Rename 'boot' column to 'unboot'
              colnames(object@summary)[colnames(object@summary) == "boot"] <- "unboot"
            }
            
            return(object)
          }
)

setMethod("pipeSubset", "ExprsPipeline",
          function(object, colBy, include){
            
            # Build an index to filter out results not matching include
            index <- object@summary[, colBy] %in% include
            
            # Subset @summary using index
            object@summary <- object@summary[index, ]
            
            # Subset @machs using index
            object@machs <- object@machs[index]
            
            return(object)
          }
)

# NOTE: IF multiple colBy terms, use product of accuracy measures
# NOTE: Increase weight of colBy elements with multiple mentions
# NOTE: IF how = 0: do not impose any threshold filter
# NOTE: IF .01 < how < 1: define threshold equal to 'how'
# NOTE: IF 1 < how < 100: define threshold equal to 'how' as a quantile
# NOTE: IF top.N = 0, include ALL @machs when building ensemble
# NOTE: IF gate != 0: impose a max filter analogous to 'how'
setMethod("pipeFilter", "ExprsPipeline",
          function(object, colBy, how = 0, gate = 0, top.N = 0){
            
            # Check if ExprsPipeline contains all colBy
            if(!all(colBy %in% colnames(object@summary))){
              
              stop("Uh oh! This ExprsPipeline does not contain all provided 'colBy' accuracy measures.")
            }
            
            # Calculate emergent top accuracy measure as product of colBy columns
            accMeasures <- apply(object@summary[colBy], MARGIN = 1, prod)
            
            # Impose threshold filter if how != 0
            if(how != 0){
              
              # Set threshold equal to 'how'
              if(0 < how & how <= 1){ threshold <- how
              
              # Set threshold equal to 'how' as a quantile
              }else if(1 < how & how <= 100){ threshold <- quantile(accMeasures, how / 100)
              
              }else if(how == "midrange"){ threshold <- (max(accMeasures) + min(accMeasures)) / 2
              
              }else if(how == "median"){ threshold <- median(accMeasures)
              
              }else if(how == "mean"){ threshold <- mean(accMeasures)
              
              }else{ stop("Uh oh! Selected 'how' not recognized!")}
              
              # Filter ExprsPipeline object
              object@summary <- object@summary[accMeasures >= threshold, ]
              object@machs <- object@machs[accMeasures >= threshold]
              
              # Filter top accuracy measure as product of colBy columns
              accMeasures <- accMeasures[accMeasures >= threshold]
            }
            
            # Impose ceiling filter if gate != 0
            if(gate != 0){
              
              # Set ceiling equal to 'gate'
              if(0 < gate & gate <= 1){ ceiling <- gate
              
              # Set ceiling equal to 'gate' as a quantile
              }else if(1 < gate & gate <= 100){ ceiling <- quantile(accMeasures, gate / 100)
              
              }else if(gate == "midrange"){ ceiling <- (max(accMeasures) + min(accMeasures)) / 2
              
              }else if(gate == "median"){ ceiling <- median(accMeasures)
              
              }else if(gate == "mean"){ ceiling <- mean(accMeasures)
              
              }else{ stop("Uh oh! Selected 'gate' not recognized!")}
              
              # Filter ExprsPipeline object
              object@summary <- object@summary[accMeasures <= ceiling, ]
              object@machs <- object@machs[accMeasures <= ceiling]
              
              # Filter top accuracy measure as product of colBy columns
              accMeasures <- accMeasures[accMeasures <= ceiling]
            }
            
            # Select top.N based on presence of a boot column
            if("boot" %in% colnames(object@summary)){
              
              # For each B boot, select 'top.N' @machs
              index <- unlist(
                
                lapply(unique(object@summary$boot),
                       function(boot){
                         
                         # Calculate total number of gridpoints for this boot
                         if(top.N > sum(object@summary$boot == boot)){
                           
                           warning("Provided 'top.N' too large for boot ", boot, ". Using all gridpoints instead.")
                           top.N <- 0
                         }
                         
                         if(top.N == 0) top.N <- sum(object@summary$boot == boot)
                         
                         # Order 'top.N' accMeasures for this boot
                         topMachs <- order(accMeasures[object@summary$boot == boot], decreasing = TRUE)[1:top.N]
                         
                         # Index by rowname for this boot
                         rownames(object@summary[object@summary$boot == boot,])[topMachs]
                       }
                )
              )
              
            }else{
              
              # Calculate total number of gridpoints for ExprsPipeline object
              if(top.N > nrow(object@summary)){
                
                warning("Provided 'top.N' too large for this ExprsPipeline object. Using all gridpoints instead.")
                top.N <- 0
              }
              
              if(top.N == 0) top.N <- nrow(object@summary)
              
              # Order 'top.N' accMeasures for entire object
              topMachs <- order(accMeasures, decreasing = TRUE)[1:top.N]
              
              # Index by rowname
              index <- rownames(object@summary)[topMachs]
            }
            
            # Filter ExprsPipeline object
            final <- rownames(object@summary) %in% index # 'final' depends on initial object@summary definition
            object@summary <- object@summary[final, ]
            object@machs <- object@machs[final]
            
            return(object)
          }
)

###########################################################
### Build ExprsEnsemble objects

# Build ExprsEnsemble using an explicit set of ExprsMachine objects
setMethod("buildEnsemble", "ExprsModel",
          function(object, ...){ # args to include additional ExprsMachine objects
            
            conjoin(object, ...)
          }
)

# Build ExprsEnsemble from an ExprsPipeline object
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
### Predict

setMethod("predict", "ExprsEnsemble",
          function(object, array, how = "probability", verbose = TRUE){
            
            if(!inherits(array, "ExprsArray")){
              
              stop("Uh oh! You can only use an ExprsEnsemble to predict on an ExprsArray object!")
            }
            
            if("ExprsMulti" %in% class(array)){
              
              stop("There does not yet exist a method for multi-class ensemble prediction!")
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
                pred[px$Case == .5] <- sample(c("Case", "Control"), sum(px$Case == .5), TRUE, prob = c(case, 1 - case))
                
                # Clean up pred
                pred <- factor(as.vector(pred), levels = c("Control", "Case"))
                names(pred) <- rownames(data)
                
              }else if(how == "majority"){
                
                # Cast a vote (1 for Case, -1 for Control)
                votes <- lapply(results, function(result) ifelse(result@pred == "Case", 1, -1))
                final <- rowSums(data.frame(votes, row.names = names(results[[1]]@pred)))
                
                # Randomly assign ties as "Case" or "Control" based on the proportion of cases
                case <- sum(array@annot$defineCase == "Case") / nrow(array@annot)
                final[final == 0] <- sample(c(1, -1), length(final[final == 0]), TRUE, prob = c(case, 1 - case))
                
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
              cat("(NOTE: This is a preview! Actual performances may differ!)\n")
              print(calcStats(final, array, aucSkip = TRUE, plotSkip = TRUE))
            }
            
            return(final)
          }
)
