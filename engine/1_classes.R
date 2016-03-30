###########################################################
### Define generic functions

setGeneric("getCases", function(object) standardGeneric("getCases"))
setGeneric("getConts", function(object) standardGeneric("getConts"))
setGeneric("getProbeSet", function(object, ...) standardGeneric("getProbeSet"))
setGeneric("getProbeSummary", function(object, ...) standardGeneric("getProbeSummary"))
setGeneric("conjoin", function(object, ...) standardGeneric("conjoin"))

###########################################################
### ExprsArray class

setClass("ExprsArray", slots = c(exprs = "matrix", annot = "data.frame", preFilter = "ANY", reductionModel = "ANY"))
setClass("ExprsBinary", contains = "ExprsArray")
setClass("ExprsMulti", contains = "ExprsArray")

setMethod("show", "ExprsArray",
          function(object){
            
            cat("##Number of classes:",
                length(unique(object@annot$defineCase)), "\n")
            
            cat("@exprs summary:",
                nrow(object@exprs), "probes by", ncol(object@exprs), "subjects\n")
            
            cat("@annot summary:",
                nrow(object@annot), "subjects by", ncol(object@annot), "annotations\n")
            
            cat("@preFilter summary:",
                unlist(lapply(object@preFilter, length)), "\n")
            
            cat("@reductionModel summary:",
                unlist(lapply(object@reductionModel, class)), "\n")
          }
)

setMethod("plot", signature(x = "ExprsArray", y = "missing"),
          function(x, i = 1, j = 2, k = 3, colors, shapes){
            
            require(lattice)
            
            # Extract components i, j, and k
            df <- data.frame(t(x@exprs))[, c(i, j, k)]
            colnames(df) <- paste0(c("i_", "j_", "k_"), colnames(df))
            
            # Plot k ~ i + j in 3D
            func <- as.formula(paste(colnames(df)[3], "~", colnames(df)[1], "+", colnames(df)[2], collapse = ""))
            if(missing(colors)) colors <- rainbow(length(unique(x@annot$defineCase)))
            if(missing(shapes)) shapes <- 19
            print(cloud(func, data = df, col = colors, pch = shapes))
          }
)

setMethod("summary", "ExprsArray",
          function(object){
            
            # For large datasets, plot only a sub-sample
            v <- as.vector(object@exprs)
            if(length(v) > 1000) v <- sample(v, 1000)
            
            # Prepare layout for multiple plots
            layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE))
            
            # Plot charts
            qqnorm(v)
            plot(density(v), main = "Density Plot")
            boxplot(v, horizontal = TRUE, main = "Box Plot", xlab = "Expression Values")
            
            # Print per-subject summary
            summary(object@exprs)
          }
)

setMethod("getProbeSet", "ExprsArray",
          function(object){
            
            return(rownames(object@exprs))
          }
)

###########################################################
### ExprsModel class

setClass("ExprsModel", slots = c(preFilter = "ANY", reductionModel = "ANY", mach = "ANY"))
setClass("ExprsMachine", contains = "ExprsModel")
setClass("ExprsModule", contains = "ExprsModel")

setMethod("show", "ExprsModel",
          function(object){
            
            cat("##Number of classes:",
                ifelse(all(class(object@mach) == "list"), length(object@mach), 2), "\n")
            
            cat("@preFilter summary:",
                unlist(lapply(object@preFilter, length)), "\n")
            
            cat("@reductionModel summary:",
                unlist(lapply(object@reductionModel, class)), "\n")
            
            cat("@mach class:", class(object@mach), "\n")
          }
)

setMethod("getProbeSet", "ExprsModel",
          function(object){
            
            return(object@preFilter[[length(object@preFilter)]])
          }
)

###########################################################
### ExprsPipeline class

setClass("ExprsPipeline", slots = c(summary = "ANY", machs = "ANY"))

setMethod("show", "ExprsPipeline",
          function(object){
            
            cat("Accuracy summary (complete summary stored in @summary slot):\n\n")
            if(nrow(object@summary) > 8){
              
              print(head(object@summary, 4))
              cat("...\n")
              print(tail(object@summary, 4))
              cat("\n")
              
            }else{
              
              print(object@summary)
              cat("\n")
            }
            
            cat("Machine summary (all machines stored in @machs slot):\n\n")
            if(length(object@machs) > 1){
              
              show(object@machs[[1]])
              cat("...\n")
              show(object@machs[[length(object@machs)]])
              cat("\n")
              
            }else{
              
              lapply(object@machs, show)
              cat("\n")
            }
          }
)

setMethod("summary", "ExprsPipeline",
          function(object){
            
            # Index which columns contain performance metrics
            index <- grepl("train.", colnames(object@summary)) | grepl("valid.", colnames(object@summary))
            
            # Calculate means and sds for performance metrics
            performance <- list("means" = apply(object@summary[, index], MARGIN = 2, mean),
                                "sds" = apply(object@summary[, index], MARGIN = 2, sd))
            
            # Tabulate parameter frequency
            parameters <- lapply(object@summary[, !index], table)
            
            # Append performance metric summary with parameter frequencies
            summary <- append(performance, parameters)
            
            return(summary)
          }
)

setMethod("getProbeSet", "ExprsPipeline",
          function(object, index){
            
            return(object@machs[[index]]@preFilter[[length(object@machs[[index]]@preFilter)]])
          }
)

setMethod("getProbeSummary", "ExprsPipeline",
          function(object){
            
            probes <- unlist(lapply(object@machs, getProbeSet))
            probes <- table(probes)[order(table(probes), decreasing = TRUE)]
            return(probes)
          }
)

###########################################################
### ExprsEnsemble class

setClass("ExprsEnsemble", slots = c(machs = "ANY"))

setMethod("show", "ExprsEnsemble",
          function(object){
            
            cat("Machine summary (all machines stored in @machs slot):\n\n")
            if(length(object@machs) > 1){
              
              show(object@machs[[1]])
              cat("...\n")
              show(object@machs[[length(object@machs)]])
              cat("\n")
              
            }else{
              
              lapply(object@machs, show)
              cat("\n")
            }
          }
)

setMethod("getProbeSet", "ExprsEnsemble",
          function(object, index){
            
            return(object@machs[[index]]@preFilter[[length(object@machs[[index]]@preFilter)]])
          }
)

setMethod("getProbeSummary", "ExprsEnsemble",
          function(object){
            
            probes <- unlist(lapply(object@machs, getProbeSet))
            probes <- table(probes)[order(table(probes), decreasing = TRUE)]
            return(probes)
          }
)

###########################################################
### ExprsPredict class

setClass("ExprsPredict", slots = c(pred = "factor", decision.values = "ANY", probability = "ANY"))

setMethod("show", "ExprsPredict",
          function(object){
            
            cat("@pred summary:", as.numeric(object@pred), "\n")
            cat("@decision.values summary:", colnames(object@decision.values), "\n")
            cat("@probability summary:", colnames(object@probability), "\n")
          }
)

###########################################################
### Misc methods

setMethod("getCases", "ExprsBinary",
          function(object){
            
            object <- new("ExprsBinary",
                          exprs = as.matrix(object@exprs[, object@annot$defineCase %in% "Case"]),
                          annot = object@annot[object@annot$defineCase %in% "Case", ],
                          preFilter = object@preFilter,
                          reductionModel = object@reductionModel)
            
            # Clean up 1-subject artifact
            colnames(object@exprs) <- rownames(object@annot)
            
            return(object)
          }
)

setMethod("getConts", "ExprsBinary",
          function(object){
            
            object <- new("ExprsBinary",
                          exprs = as.matrix(object@exprs[, object@annot$defineCase %in% "Control"]),
                          annot = object@annot[object@annot$defineCase %in% "Control", ],
                          preFilter = object@preFilter,
                          reductionModel = object@reductionModel)
            
            # Clean up 1-subject artifact
            colnames(object@exprs) <- rownames(object@annot)
            
            return(object)
          }
)

###########################################################
### Conjoin

# Combines multiple ExprsArray objects into single ExprsArray object
# NOTE: Function only works on ExprsArray objects that have not undergone feature selection
# NOTE: Any annotations missing in an @annot replaced with NA values
setMethod("conjoin", "ExprsArray",
          function(object, ...){ # args to include additional ExprsArray objects
            
            require(plyr)
            
            # Prepare list of ExprsArray objects
            args <- list(...)
            index <- unlist(lapply(args, function(arg) inherits(arg, "ExprsArray")))
            args <- append(list(object), args[index])
            
            if(!length(args) > 1 | any(!sapply(args, function(e) identical(class(e), class(object))))){
              
              stop("User must provide additional ExprsArray objects of the same class!")
            }
            
            if(any(!sapply(args, function(e) is.null(e@preFilter) | is.null(e@reductionModel)))){
              
              stop("This function is not equipped to handle objects that have undergone feature selection!")
            }
            
            # Prepare single matrix for @exprs and @annot each
            exprs <- as.matrix(do.call(cbind, lapply(args, function(a) a@exprs)))
            annot <- do.call(rbind.fill, lapply(args, function(a) a@annot))
            rownames(annot) <- unlist(lapply(args, function(a) rownames(a@annot)))
            
            # Return single ExprsArray object
            new(class(object), exprs = exprs, annot = annot, preFilter = NULL, reductionModel = NULL)
          }
)

setMethod("conjoin", "ExprsModel",
          function(object, ...){ # args to include additional ExprsModel objects
            
            # Prepare list of ExprsModel objects
            args <- list(...)
            index <- unlist(lapply(args, function(arg) inherits(arg, "ExprsModel")))
            machs <- append(object, args[index])
            
            # Return single ExprsEnsemble object
            new("ExprsEnsemble", machs = machs)
          }
)

setMethod("conjoin", "ExprsPipeline",
          function(object, ...){ # args to include additional ExprsPipeline objects
            
            require(plyr)
            
            # Prepare list of ExprsPipeline objects
            args <- list(...)
            index <- unlist(lapply(args, function(arg) class(arg) == "ExprsPipeline"))
            args.summary <- append(list(object@summary), lapply(args[index], function(pl) pl@summary))
            args.machs <- append(list(object@machs), lapply(args[index], function(pl) pl@machs))
            
            # Initialize the conjoin boot counter
            b <- 1
            
            # Apply conjoin boot counter to each ExprsPipeline object
            pls <- lapply(args.summary,
                          function(pl){
                            
                            # Initialize the conjoin boot container
                            pl <- cbind("join" = 0, pl)
                            
                            if(!"boot" %in% colnames(pl)){
                              
                              # Add conjoin boot counter
                              pl$join <- b
                              b <<- b + 1
                              
                            }else{
                              
                              # For each boot in $boot
                              for(i in 1:length(unique(pl$boot))){
                                
                                # Change each unique boot to conjoin boot counter
                                pl$join[pl$boot == i] <- b
                                b <<- b + 1
                              }
                              
                              # Rename $boot to $unboot
                              colnames(pl)[colnames(pl) == "boot"] <- "unboot"
                            }
                            
                            # Rename $join to $boot
                            colnames(pl)[colnames(pl) == "join"] <- "boot"
                            
                            return(pl)
                          }
            )
            
            # Return single ExprsPipeline object
            new("ExprsPipeline", summary = do.call(rbind.fill, pls), machs = unlist(args.machs))
          }
)

setMethod("conjoin", "ExprsEnsemble",
          function(object, ...){ # args to include additional ExprsEnsemble objects
            
            # Prepare list of ExprsEnsemble objects
            args <- list(...)
            index <- unlist(lapply(args, function(arg) class(arg) == "ExprsEnsemble"))
            machs <- unlist(append(list(object@machs), lapply(args[index], function(pl) pl@machs)))
            
            # Return single ExprsEnsemble object
            new("ExprsEnsemble", machs = machs)
          }
)
