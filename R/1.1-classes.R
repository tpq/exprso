###########################################################
### Define generic functions

#' Retrieve \code{exprso} Feature Set
#'
#' See the respective S4 class for method details.
#'
#' @param object An \code{ExprsArray}, \code{ExprsModel}, \code{ExprsPipeline},
#'  \code{ExprsEnsemble}, or \code{ExprsPredict} object.
#' @param ... See \code{\link{ExprsPipeline-class}} or \code{\link{ExprsEnsemble-class}}.
#'
#' @seealso
#' * \code{\link{ExprsArray-class}}
#' * \code{\link{ExprsModel-class}}
#' * \code{\link{ExprsPipeline-class}}
#' * \code{\link{ExprsEnsemble-class}}
#' * \code{\link{ExprsPredict-class}}
#' @export
setGeneric("getProbeSet",
           function(object, ...) standardGeneric("getProbeSet")
)

###########################################################
### ExprsArray class

#' An S4 class to store feature and annotation data
#'
#' @slot exprs A matrix. Stores the expression data.
#' @slot annot A data.frame. Stores the annotation data.
#' @slot preFilter Typically a list. Stores feature selection history.
#' @slot reductionModel Typically a list. Stores dimension reduction history.
#'
#' @seealso
#' * \code{\link{ExprsArray-class}}
#' * \code{\link{ExprsModel-class}}
#' * \code{\link{ExprsPipeline-class}}
#' * \code{\link{ExprsEnsemble-class}}
#' * \code{\link{ExprsPredict-class}}
#' @export
setClass("ExprsArray",
         slots = c(
           exprs = "matrix",
           annot = "data.frame",
           preFilter = "ANY",
           reductionModel = "ANY"
         )
)

#' @describeIn ExprsArray-class An \code{ExprsArray} subclass for data with binary class labels.
#' @export
setClass("ExprsBinary", contains = "ExprsArray")

#' @describeIn ExprsArray-class An \code{ExprsArray} subclass for data with many class labels.
#' @export
setClass("ExprsMulti", contains = "ExprsArray")

#' @describeIn ExprsArray-class Method to show \code{ExprsArray} object.
#'
#' @param object,x An object of class \code{ExprsArray}.
#' @export
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

#' @describeIn ExprsArray-class Method to subset \code{ExprsArray} object.
#'
# #' @param x An object of class \code{ExprsArray}.
#' @param i,j,drop Subsets via \code{object@annot[i, j, drop]}.
#' @export
setMethod('[', signature(x = "ExprsArray"),
          function(x, i, j, drop){

            if(!missing(j)){

              return(x@annot[i, j, drop])

            }else{

              x@annot <- x@annot[i, j, drop]
              x@exprs <- x@exprs[, rownames(x@annot), drop]
              return(x)
            }
          }
)

#' @describeIn ExprsArray-class Method to subset \code{ExprsArray} object.
#'
# #' @param x An object of class \code{ExprsArray}.
#' @param name Subsets via \code{object@annot[, name]}.
#' @export
setMethod('$', signature(x = "ExprsArray"),
          function(x, name){

            return(x@annot[, name])
          }
)

#' @describeIn ExprsArray-class Method to plot three dimensions of the expression data.
#'
#' @param i,j,k A numeric scalar. Indexes the first, second, and third dimensions to plot.
#' @param colors A character vector. Optional. Manually assign a color to each subject point.
#' @param shapes A numeric vector. Optional. Manually assign a shape to each subject point.
#'
#' @import lattice
#' @export
setMethod("plot", signature(x = "ExprsArray", y = "missing"),
          function(x, i = 1, j = 2, k = 3, colors, shapes){

            # Extract components i, j, and k
            df <- data.frame(t(x@exprs))[, c(i, j, k)]
            colnames(df) <- paste0(c("i_", "j_", "k_"), colnames(df))

            # Plot k ~ i + j in 3D
            func <- as.formula(paste(colnames(df)[3], "~",
                                     colnames(df)[1], "+",
                                     colnames(df)[2],
                                     collapse = ""))
            if(missing(colors)) colors <- rainbow(length(unique(x@annot$defineCase)))
            if(missing(shapes)) shapes <- 19
            print(lattice::cloud(func, data = df, col = colors, pch = shapes))
          }
)

#' @describeIn ExprsArray-class Method to plot summary graphs for a sub-sample of expression data.
#' @export
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

            layout(matrix(c(1), 1, 1, byrow = TRUE))

            # Print per-subject summary
            summary(object@exprs)
          }
)

#' @describeIn ExprsArray-class Method to return features within an \code{ExprsArray} object.
#' @export
setMethod("getProbeSet", "ExprsArray",
          function(object){

            return(rownames(object@exprs))
          }
)

###########################################################
### ExprsModel class

#' An S4 class to store the classification model
#'
#' @slot preFilter Typically a list. Stores feature selection history.
#' @slot reductionModel Typically a list. Stores dimension reduction history.
#' @slot mach Typically an S4 class. Stores the classification model.
#'
#' @seealso
#' * \code{\link{ExprsArray-class}}
#' * \code{\link{ExprsModel-class}}
#' * \code{\link{ExprsPipeline-class}}
#' * \code{\link{ExprsEnsemble-class}}
#' * \code{\link{ExprsPredict-class}}
#' @export
setClass("ExprsModel",
         slots = c(
           preFilter = "ANY",
           reductionModel = "ANY",
           mach = "ANY"
         )
)

#' @describeIn ExprsModel-class An \code{ExprsModel} subclass for dichotomous classifiers.
#' @export
setClass("ExprsMachine", contains = "ExprsModel")

#' @describeIn ExprsModel-class An \code{ExprsModel} subclass for multi-class classifiers.
#' @export
setClass("ExprsModule", contains = "ExprsModel")

#' @describeIn ExprsModel-class Method to show \code{ExprsModel} object.
#'
#' @param object An object of class \code{ExprsModel}.
#' @export
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

#' @describeIn ExprsModel-class Method to return features within an \code{ExprsModel} object.
#' @export
setMethod("getProbeSet", "ExprsModel",
          function(object){

            return(object@preFilter[[length(object@preFilter)]])
          }
)

###########################################################
### ExprsPipeline class

#' An S4 class to store models built during high-throughput learning
#'
#' @slot summary Typically a data.frame. Stores the parameters and performances for classification models.
#' @slot machs Typically a list. Stores the classification models referenced in \code{summary} slot.
#'
#' @seealso
#' * \code{\link{ExprsArray-class}}
#' * \code{\link{ExprsModel-class}}
#' * \code{\link{ExprsPipeline-class}}
#' * \code{\link{ExprsEnsemble-class}}
#' * \code{\link{ExprsPredict-class}}
#' @export
setClass("ExprsPipeline", slots = c(summary = "ANY", machs = "ANY"))

#' @describeIn ExprsPipeline-class Method to show \code{ExprsPipeline} object.
#'
#' @param object An object of class \code{ExprsPipeline}.
#' @export
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

#' @describeIn ExprsPipeline-class Method to subset \code{ExprsPipeline} object.
#'
# #' @param x An object of class \code{ExprsPipeline}.
#' @param i,j,drop Subsets via \code{object@summary[i, j, drop]}.
#' @export
setMethod('[', signature(x = "ExprsPipeline"),
          function(x, i, j, drop){

            if(!missing(j)){

              return(x@summary[i, j, drop])

            }else{

              index <- which(rownames(x@summary) %in% rownames(x@summary[i, j, drop]))
              x@summary <- x@summary[index, j, drop]
              x@machs <- x@machs[index]
              return(x)
            }
          }
)

#' @describeIn ExprsPipeline-class Method to subset \code{ExprsPipeline} object.
#'
# #' @param x An object of class \code{ExprsPipeline}.
#' @param name Subsets via \code{object@summary[, name]}.
#' @export
setMethod('$', signature(x = "ExprsPipeline"),
          function(x, name){

            return(x@summary[, name])
          }
)

#' @describeIn ExprsPipeline-class Method to summarize \code{ExprsPipeline} parameters and performances.
#' @export
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

#' @describeIn ExprsPipeline-class Method to return features within an \code{ExprsPredict} model.
#'
#' @param index A numeric scalar. Indicates the i-th model from which to retrieve the features. If missing,
#'  \code{getProbeSet} will tabulate features across all models.
#' @export
setMethod("getProbeSet", "ExprsPipeline",
          function(object, index){

            if(!missing(index)){

              return(object@machs[[index]]@preFilter[[length(object@machs[[index]]@preFilter)]])

            }else{

              probes <- unlist(lapply(object@machs, getProbeSet))
              probes <- table(probes)[order(table(probes), decreasing = TRUE)]
              return(probes)
            }
          }
)

###########################################################
### ExprsEnsemble class

#' An S4 class to store multiple classification models
#'
#' @inheritParams ExprsPipeline-class
#'
#' @seealso
#' * \code{\link{ExprsArray-class}}
#' * \code{\link{ExprsModel-class}}
#' * \code{\link{ExprsPipeline-class}}
#' * \code{\link{ExprsEnsemble-class}}
#' * \code{\link{ExprsPredict-class}}
#' @export
setClass("ExprsEnsemble",
         slots = c(
           machs = "ANY"
         )
)

#' @describeIn ExprsEnsemble-class Method to show \code{ExprsEnsemble} object.
#'
#' @param object An object of class \code{ExprsEnsemble}.
#' @export
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

#' @describeIn ExprsEnsemble-class Method to return features within an \code{ExprsEnsemble} model.
#'
#' @inheritParams ExprsPipeline-class
#' @export
setMethod("getProbeSet", "ExprsPipeline",
          function(object, index){

            if(!missing(index)){

              return(object@machs[[index]]@preFilter[[length(object@machs[[index]]@preFilter)]])

            }else{

              probes <- unlist(lapply(object@machs, getProbeSet))
              probes <- table(probes)[order(table(probes), decreasing = TRUE)]
              return(probes)
            }
          }
)

###########################################################
### ExprsPredict class

#' An S4 class to store class predictions
#'
#' @slot pred A factor. Stores class predictions as an unambiguous class assignment.
#' @slot decision.values Typically a matrix. Stores class predictions as a decision value.
#' @slot probability Typically a matrix. Stores class predictions as a probability.
#'
#' @seealso
#' * \code{\link{ExprsArray-class}}
#' * \code{\link{ExprsModel-class}}
#' * \code{\link{ExprsPipeline-class}}
#' * \code{\link{ExprsEnsemble-class}}
#' * \code{\link{ExprsPredict-class}}
#' @export
setClass("ExprsPredict",
         slots = c(
           pred = "factor",
           decision.values = "ANY",
           probability = "ANY"
         )
)

#' @describeIn ExprsPredict-class Method to show \code{ExprsPredict} object.
#'
#' @param object An object of class \code{ExprsPredict}.
#' @export
setMethod("show", "ExprsPredict",
          function(object){

            cat("@pred summary:", as.numeric(object@pred), "\n")
            cat("@decision.values summary:", colnames(object@decision.values), "\n")
            cat("@probability summary:", colnames(object@probability), "\n")
          }
)
