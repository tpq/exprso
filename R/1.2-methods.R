###########################################################
### Define generic functions

#' Retrieve Feature Set
#'
#' See the respective S4 class for method details.
#'
#' @param object An \code{ExprsArray}, \code{ExprsModel}, \code{ExprsPipeline},
#'  or \code{ExprsEnsemble} object.
#' @param ... See \code{\link{ExprsPipeline-class}} or
#'  \code{\link{ExprsEnsemble-class}}.
#'
#' @seealso
#' \code{\link{ExprsArray-class}}\cr
#' \code{\link{ExprsModel-class}}\cr
#' \code{\link{ExprsPipeline-class}}\cr
#' \code{\link{ExprsEnsemble-class}}\cr
#' \code{\link{ExprsPredict-class}}
#' @export
setGeneric("getFeatures",
           function(object, ...) standardGeneric("getFeatures")
)

###########################################################
### ExprsArray class

#' @describeIn ExprsArray Method to show \code{ExprsArray} object.
#'
#' @param object,x An object of class \code{ExprsArray}.
#'
#' @export
setMethod("show", "ExprsArray",
          function(object){

            cat("##Number of classes:",
                length(unique(object@annot$defineCase)), "\n")

            cat("@exprs summary:",
                nrow(object@exprs), "features by", ncol(object@exprs), "subjects\n")

            cat("@annot summary:",
                nrow(object@annot), "subjects by", ncol(object@annot), "annotations\n")

            cat("@preFilter summary:",
                unlist(lapply(object@preFilter, length)), "\n")

            cat("@reductionModel summary:",
                unlist(lapply(object@reductionModel, class)), "\n")
          }
)

#' @describeIn ExprsArray Method to subset \code{ExprsArray} object.
#'
#' @param i,j Subsets entire \code{ExprsArray} object via
#'  \code{object@annot[i, j]}. Returns \code{object@annot[, j]} if
#'  argument \code{i} is missing.
#'
#' @aliases [,ExprsArray-method
#' @docType methods
#' @export
setMethod('[', signature(x = "ExprsArray", i = "ANY", j = "ANY"),
          function(x, i, j){

            if(!missing(j)){

              return(x@annot[i, j])

            }else{

              x@annot <- as.data.frame(x@annot[i, j, drop = FALSE])
              x@exprs <- x@exprs[, i, drop = FALSE]
              colnames(x@exprs) <- rownames(x@annot)
              return(x)
            }
          }
)

#' @describeIn ExprsArray Method to subset \code{ExprsArray} object.
#'
#' @param name Returns \code{object@annot[, name]}.
#'
#' @export
setMethod('$', signature(x = "ExprsArray"),
          function(x, name){

            return(x@annot[, name])
          }
)

#' @describeIn ExprsArray Method to subset \code{ExprsArray} object.
#'
#' @param subset Subsets entire \code{ExprsArray} object via
#'  \code{object@annot[subset, ]}. Can be used to rearrange feature order.
#' @param select Subsets entire \code{ExprsArray} object via
#'  \code{object@annot[, select]}. Can be used to rearrange subject order.
#'
#' @export
setMethod("subset", signature(x = "ExprsArray"),
          function(x, subset, select){

            if(missing(subset)) subset <- rownames(x@annot)
            if(missing(select)) select <- colnames(x@annot)

            x@annot <- x@annot[subset, select, drop = FALSE]
            x@exprs <- x@exprs[, subset, drop = FALSE]

            return(x)
          }
)

#' @describeIn ExprsArray Method to quickly plot two or three dimensions of data.
#'
#' @param y Leave missing. Argument exists because of \code{\link{plot}} generic definition.
#' @param a,b,c A numeric scalar. Indexes the first, second, and third dimensions to plot.
#'  Set \code{c = 0} to plot two dimensions.
#' @param ... Additional arguments passed to\code{plot} or \code{lattice::cloud}.
#'
#' @importFrom stats as.formula
#' @importFrom grDevices rainbow
#' @importFrom lattice cloud
#' @export
setMethod("plot", signature(x = "ExprsArray", y = "missing"),
          function(x, y, a = 1, b = 2, c = 3, ...){

            args <- getArgs(...)

            args <- defaultArg("col", grDevices::rainbow(length(unique(x$defineCase)))
                               [as.numeric(factor(x$defineCase))], args)
            args <- defaultArg("pch", 19, args)

            if(c > 0){

              # Extract components a, b, and c
              df <- data.frame(t(x@exprs))[, c(a, b, c)]
              colnames(df) <- paste0(c("a_", "b_", "c_"), colnames(df))

              # Plot c ~ a + b in 3D
              func <- stats::as.formula(paste(colnames(df)[3], "~",
                                              colnames(df)[1], "+",
                                              colnames(df)[2],
                                              collapse = ""))

              args <- append(args, list("x" = func, "data" = df))
              do.call("cloud", args)

            }else{

              # Extract components a and b
              df <- data.frame(t(x@exprs))[, c(a, b)]
              colnames(df) <- paste0(c("a_", "b_"), colnames(df))

              args <- defaultArg("xlab", colnames(df)[1], args)
              args <- defaultArg("ylab", colnames(df)[2], args)

              args <- append(args, list("x" = df[, 1], "y" = df[, 2]))
              do.call("plot", args)
            }
          }
)

#' @describeIn ExprsArray Method to plot summary graphs for a sub-sample of feature data.
#'
#' @importFrom stats qqnorm density
#' @importFrom graphics layout boxplot
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
            boxplot(v, horizontal = TRUE, main = "Box Plot", xlab = "Feature Values")

            layout(matrix(c(1), 1, 1, byrow = TRUE))

            # Print per-subject summary
            summary(object@exprs)
          }
)

#' @describeIn ExprsArray Method to return features within an \code{ExprsArray} object.
#'
#' @export
setMethod("getFeatures", "ExprsArray",
          function(object){

            return(rownames(object@exprs))
          }
)

###########################################################
### ExprsModel class

#' @describeIn ExprsModel Method to show \code{ExprsModel} object.
#'
#' @param object An object of class \code{ExprsModel}.
#'
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

#' @describeIn ExprsModel Method to return features within an \code{ExprsModel} object.
#'
#' @export
setMethod("getFeatures", "ExprsModel",
          function(object){

            return(object@preFilter[[length(object@preFilter)]])
          }
)

###########################################################
### ExprsPipeline class

#' @describeIn ExprsPipeline Method to show \code{ExprsPipeline} object.
#'
#' @param object,x An object of class \code{ExprsPipeline}.
#'
#' @importFrom utils head tail
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

#' @describeIn ExprsPipeline Method to subset \code{ExprsPipeline} object.
#'
#' @param i,j Subsets entire \code{ExprsPipeline} object via
#'  \code{object@summary[i, j]}. Returns \code{object@summary[, j]} if
#'  argument \code{i} is missing.
#'
#' @aliases [,ExprsPipeline-method
#' @docType methods
#' @export
setMethod('[', signature(x = "ExprsPipeline", i = "ANY", j = "ANY"),
          function(x, i, j){

            if(!missing(j)){

              return(x@summary[i, j])

            }else{

              index <- which(rownames(x@summary) %in% rownames(x@summary[i, j, drop = FALSE]))
              x@summary <- x@summary[index, j, drop = FALSE]
              x@machs <- x@machs[index]
              return(x)
            }
          }
)

#' @describeIn ExprsPipeline Method to subset \code{ExprsPipeline} object.
#'
#' @param name Returns \code{object@summary[, name]}.
#'
#' @export
setMethod('$', signature(x = "ExprsPipeline"),
          function(x, name){

            return(x@summary[, name])
          }
)

#' @describeIn ExprsPipeline Method to subset \code{ExprsPipeline} object.
#'
#' @param subset Subsets entire \code{ExprsPipeline} object via
#'  \code{object@summary[subset, ]}. Can be used to rearrange summary table.
#' @param select Subsets entire \code{ExprsPipeline} object via
#'  \code{object@summary[, select]}. Can be used to rearrange summary table.
#'
#' @export
setMethod("subset", signature(x = "ExprsPipeline"),
          function(x, subset, select){

            if(missing(subset)) subset <- rownames(x@summary)
            if(missing(select)) select <- colnames(x@summary)

            x@summary <- x@summary[subset, select, drop = FALSE]
            x@machs <- x@machs[subset]

            return(x)
          }
)

#' @describeIn ExprsPipeline Method to summarize \code{ExprsPipeline} classification results.
#'
#' @importFrom stats sd
#' @export
setMethod("summary", "ExprsPipeline",
          function(object){

            # Index which columns contain performance metrics
            index <- grepl("train.", colnames(object@summary)) | grepl("valid.", colnames(object@summary))

            # Calculate means and sds for performance metrics
            performance <- list("means" = apply(object@summary[, index], MARGIN = 2, mean),
                                "sds" = apply(object@summary[, index], MARGIN = 2, sd))

            # Tabulate parameter frequency
            parameters <- lapply(object@summary[, !index], function(column) table(as.character(column)))

            # Append performance metric summary with parameter frequencies
            summary <- append(performance, parameters)

            return(summary)
          }
)

#' @describeIn ExprsPipeline Method to return features within an \code{ExprsPredict} model.
#'
#' @param index A numeric scalar. The i-th model from which to retrieve features.
#'  If missing, \code{getFeatures} will tabulate features across all models.
#'
#' @export
setMethod("getFeatures", "ExprsPipeline",
          function(object, index){

            if(!missing(index)){

              return(object@machs[[index]]@preFilter[[length(object@machs[[index]]@preFilter)]])

            }else{

              features <- unlist(lapply(object@machs, getFeatures))
              features <- table(features)[order(table(features), decreasing = TRUE)]
              return(features)
            }
          }
)

###########################################################
### ExprsEnsemble class

#' @describeIn ExprsEnsemble Method to show \code{ExprsEnsemble} object.
#'
#' @param index A numeric scalar. The i-th model from which to retrieve features.
#'  If missing, \code{getFeatures} will tabulate features across all models.
#'
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

#' @describeIn ExprsEnsemble Method to return features within an \code{ExprsEnsemble} model.
#'
#' @inheritParams getFeatures
#'
#' @export
setMethod("getFeatures", "ExprsEnsemble",
          function(object, index){

            if(!missing(index)){

              return(object@machs[[index]]@preFilter[[length(object@machs[[index]]@preFilter)]])

            }else{

              features <- unlist(lapply(object@machs, getFeatures))
              features <- table(features)[order(table(features), decreasing = TRUE)]
              return(features)
            }
          }
)

###########################################################
### ExprsPredict class

#' @describeIn ExprsPredict Method to show \code{ExprsPredict} object.
#'
#' @param object An object of class \code{ExprsPredict}.
#'
#' @export
setMethod("show", "ExprsPredict",
          function(object){

            cat("@pred summary:", as.numeric(object@pred), "\n")
            cat("@decision.values summary:", colnames(object@decision.values), "\n")
            cat("@probability summary:", colnames(object@probability), "\n")
          }
)
