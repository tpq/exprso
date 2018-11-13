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

            if(class(object) == "ExprsBinary" | class(object) == "ExprsMulti"){
              cat("##Number of classes:",
                  length(unique(object@annot$defineCase)), "\n")
            }else{
              cat("##Average outcome: ",
                  round(mean(object@annot$defineCase), 2), " [",
                  round(min(object@annot$defineCase), 2), "-",
                  round(max(object@annot$defineCase), 2), "]\n", sep = "")
            }

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

#' @describeIn ExprsArray Method to plot two or three dimensions of data.
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

            classCheck(x, c("ExprsBinary", "ExprsMulti"),
                       "This plot method only works for classification tasks.")

            args <- getArgs(...)
            args <- defaultArg("pch", 19, args)
            args <- defaultArg("col", grDevices::rainbow(length(unique(x$defineCase)))
                               [as.numeric(factor(x$defineCase))], args)

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

            if(class(object) == "ExprsMachine" | class(object) == "ExprsModule"){
              cat("##Number of classes:",
                  ifelse(all(class(object@mach) == "list"), length(object@mach), 2), "\n")
            }else{
              cat("##Continuous outcome model\n")
            }

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

#' @describeIn ExprsPipeline Method to summarize \code{ExprsPipeline} results.
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
#' @param index A numeric scalar. The i-th model from which to retrieve features or weights.
#'  If missing, function will tabulate features or weights across all models.
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
#' @param index A numeric scalar. The i-th model from which to retrieve features or weights.
#'  If missing, function will tabulate features or weights across all models.
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

            cat("@pred summary:", table(as.numeric(object@pred)), "\n")
            cat("@decision.values summary:", colnames(object@decision.values), "\n")
            cat("@probability summary:", colnames(object@probability), "\n")
            if(length(levels(object@pred)) == 2){
              cat("@actual summary:", table(as.numeric(
                factor(object@actual, levels = levels(object@pred)))), "\n")
            }else{
              cat("@actual summary:", table(as.numeric(
                object@actual)), "\n")
            }
          }
)

#' @describeIn RegrsPredict Method to show \code{RegrsPredict} object.
#'
#' @param object An object of class \code{RegrsPredict}.
#'
#' @export
setMethod("show", "RegrsPredict",
          function(object){

            cat("@pred summary: ",
                round(mean(object@pred), 2), " [",
                round(min(object@pred), 2), "-",
                round(max(object@pred), 2), "]\n", sep = "")
            cat("@actual summary: ",
                round(mean(object@actual), 2), " [",
                round(min(object@actual), 2), "-",
                round(max(object@actual), 2), "]\n", sep = "")
          }
)

###########################################################
### Tidy wrappers

#' Extract Training Set
#'
#' This function extracts the training set from the result of a
#'  \code{split} method call such as \code{splitSample} or \code{splitStratify}.
#'
#' @param splitSets A two-item list. The result of a \code{split} method call.
#' @return An \code{ExprsArray} object.
#' @export
trainingSet <- function(splitSets){

  if(class(splitSets) == "list" & length(splitSets) == 2){
    return(splitSets[["array.train"]])
  }else{
    stop("Uh oh! Cannot extract a training set from this object.")
  }
}

#' Extract Validation Set
#'
#' This function extracts the validation set from the result of a
#'  \code{split} method call such as \code{splitSample} or \code{splitStratify}.
#'
#' @inheritParams trainingSet
#' @return An \code{ExprsArray} object.
#' @export
validationSet <- function(splitSets){

  if(class(splitSets) == "list" & length(splitSets) == 2){
    return(splitSets[["array.valid"]])
  }else{
    stop("Uh oh! Cannot extract a test set from this object.")
  }
}

#' @describeIn validationSet A variant of \code{validationSet}.
#' @export
testSet <- function(splitSets){

  validationSet(splitSets)
}

#' Tidy Subset Wrapper
#'
#' \code{modSubset} function provides a tidy wrapper for the \code{ExprsArray}
#'  \code{subset} method. \code{pipeSubset} provides a tidy wrapper for the
#'  \code{ExprsPipeline} \code{subset} method.
#'
#' @inheritParams arrayExprs
#' @param object An \code{ExprsArray} or \code{ExprsPipeline} object to subset.
#' @param include A character vector. Specifies which annotations in \code{colBy}
#'  to include in the subset.
#' @return An \code{ExprsArray} or \code{ExprsPipeline} object.
#' @export
modSubset <- function(object, colBy, include){

  if(inherits(object, "ExprsArray") | class(object) == "ExprsPipeline"){
    subset(object, subset = object[, colBy] %in% include)
  }else{
    stop("Uh oh! You can only use modSubset on an ExprsArray or ExprsPipeline object!")
  }
}

#' @describeIn modSubset A variant of \code{modSubset}.
#' @export
pipeSubset <- function(object, colBy, include){

  modSubset(object, colBy, include)
}

###########################################################
### getWeights LASSO methods

#' Retrieve LASSO Weights
#'
#' See the respective S4 class for method details.
#'
#' @param object An \code{ExprsModel}, \code{ExprsPipeline},
#'  or \code{ExprsEnsemble} object.
#' @param ... Arguments passed to \code{glmnet::coef.cv.glmnet}.
#'
#' @export
setGeneric("getWeights",
           function(object, ...) standardGeneric("getWeights")
)

#' @describeIn ExprsModel Method to return LASSO weights.
#'
#' @param ... For \code{getWeights}, optional arguments passed to
#'  \code{glmnet::coef.cv.glmnet}.
#'
#' @export
setMethod("getWeights", "ExprsModel",
          function(object, ...){

            if("cv.glmnet" %in% class(object@mach)){

              df <- t(as.matrix(glmnet::coef.cv.glmnet(object@mach, ...)))
              df <- data.frame(df)
              colnames(df)[1] <- "Intercept"
              return(df)

            }else if("randomForest" %in% class(object@mach)){

              catch <- randomForest::importance(object@mach, ...)
              how <- colnames(catch)
              df <- t(catch[,1])
              df <- data.frame("importance" = how, df)
              return(df)

            }else{

              stop("This method only works for 'cv.glmnet' or 'randomForest' objects.")
            }
          }
)

#' @describeIn ExprsPipeline Method to return LASSO weights.
#'
#' @param ... For \code{getWeights}, optional arguments passed to
#'  \code{glmnet::coef.cv.glmnet}.
#'
#' @export
setMethod("getWeights", "ExprsPipeline",
          function(object, index, ...){

            if(!missing(index)){

              getWeights(object@machs[[index]])

            }else{

              weights <- lapply(object@machs, getWeights)
              weights.df <- do.call(plyr::rbind.fill, weights)
              final <- cbind(object@summary, weights.df)
              return(final)
            }
          }
)

#' @describeIn ExprsEnsemble Method to return LASSO weights.
#'
#' @param ... For \code{getWeights}, optional arguments passed to
#'  \code{glmnet::coef.cv.glmnet}.
#'
#' @export
setMethod("getWeights", "ExprsEnsemble",
          function(object, index, ...){

            if(!missing(index)){

              getWeights(object@machs[[index]])

            }else{

              weights <- lapply(object@machs, getWeights)
              weights.df <- do.call(plyr::rbind.fill, weights)
              return(weights.df)
            }
          }
)
