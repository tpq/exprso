#' Workhorse for build Methods
#'
#' Used as a back-end wrapper for creating new build methods.
#'
#' @inheritParams fs.
#' @param object An \code{ExprsArray} object. The training set.
#' @return Returns an \code{ExprsModel} object.
#' @export
build. <- function(object, top, uniqueFx, ...){

  # Convert top input to explicit feature reference
  if(class(top) == "numeric"){
    if(length(top) == 1){
      if(top > nrow(object@exprs)) top <- 0
      if(top == 0){
        top <- rownames(object@exprs) # keep for uniqueFx
        if(class(object) != "ExprsMulti") data <- t(object@exprs)
      }else{
        top <- rownames(object@exprs)[1:top]
        if(class(object) != "ExprsMulti") data <- t(object@exprs[top, , drop = FALSE])
      }
    }else{
      top <- rownames(object@exprs)[top] # when top is numeric vector
      if(class(object) != "ExprsMulti") data <- t(object@exprs[top, , drop = FALSE])
    }
  }else{
    # top <- top # feature names already specified
    if(class(object) != "ExprsMulti") data <- t(object@exprs[top, , drop = FALSE])
  }

  # Carry through and append fs history as stored in the ExprsArray object
  # NOTE: length(ExprsMachine@preFilter) > length(ExprsArray@preFilter)
  if(class(object) == "ExprsMulti"){
    args <- getArgs(...)
    args <- append(args, list("object" = object, "top" = top, "method" = "build.", "uniqueFx" = uniqueFx))
    machs <- do.call("doMulti", args)
    m <- new("ExprsModule",
             preFilter = append(object@preFilter, list(top)),
             reductionModel = append(object@reductionModel, list(NA)),
             mach = machs)
  }else if(class(object) == "ExprsBinary"){
    labels <- factor(object@annot$defineCase, levels = c("Control", "Case"))
    model <- do.call("uniqueFx", list(data, labels, ...))
    m <- new("ExprsMachine",
             preFilter = append(object@preFilter, list(top)),
             reductionModel = append(object@reductionModel, list(NA)),
             mach = model)
  }else if(class(object) == "RegrsArray"){
    labels <- object@annot$defineCase
    model <- do.call("uniqueFx", list(data, labels, ...))
    m <- new("RegrsModel",
             preFilter = append(object@preFilter, list(top)),
             reductionModel = append(object@reductionModel, list(NA)),
             mach = model)
  }else{ stop("Uh oh! DEBUG ERROR: 003")}

  return(m)
}

#' Perform Multiple "1 vs. all" Tasks
#'
#' A function to execute multiple "1 vs. all" binary tasks.
#'
#' \code{doMulti} runs once for each factor level in the
#'  "defineCase" column. If a training set is missing any
#'  one of the factor levels (e.g., owing to random cuts during
#'  cross-validation), the \code{ExprsModule} component that
#'  would refer to that class label gets replaced with an NA
#'  placeholder. Note that this NA placeholder will prevent a
#'  classifier from possibly predicting the NA class (i.e., a
#'  classifier can only make predictions about class
#'  labels that it "knows"). However, these "unknown" classes
#'  still impact metrics of classifier performance.
#'  Otherwise, see \code{\link{exprso-predict}}.
#'
#' @inheritParams build.
#' @param method A character string. The method to apply.
#' @return A list of the results from \code{method}.
#' @export
doMulti <- function(object, top = 0, method, ...){

  classCheck(object, "ExprsMulti",
             "This method only works for multi-class classification tasks.")

  # Perform N binary tasks
  args <- getArgs(...)
  multi <- vector("list", length(levels(object@annot$defineCase)))
  for(i in 1:length(levels(object@annot$defineCase))){

    # If the i-th ExprsMachine would not have any representative cases
    if(all(as.numeric(object@annot$defineCase) != i)){

      cat("Missing class ", i, ". Using a NA placeholder instead.\n", sep = "")
      multi[[i]] <- NA

    }else{

      # Turn the ExprsMulti object into the i-th ExprsBinary object
      temp <- object
      temp@annot$defineCase <- ifelse(as.numeric(temp$defineCase) == i, "Case", "Control")
      class(temp) <- "ExprsBinary"

      # Perform the binary task
      cat("Performing a one-vs-all binary task with class", i, "set as \"Case\".\n")
      args.i <- append(list("object" = temp, "top" = top), args)
      multi[[i]] <- do.call(method, args.i)
    }
  }

  return(multi)
}

#' Build Naive Bayes Model
#'
#' \code{buildNB} builds a model using the \code{naiveBayes} function
#'  from the \code{e1071} package.
#'
#' @inheritParams build.
#' @return Returns an \code{ExprsModel} object.
#' @export
buildNB <- function(object, top = 0, ...){ # args to naiveBayes

  classCheck(object, c("ExprsBinary", "ExprsMulti"),
             "This build method only works for classification tasks.")

  build.(object, top,
         uniqueFx = function(data, labels, ...){

           # Perform naiveBayes via ~ method
           args <- getArgs(...)
           df <- data.frame(data, "defineCase" = labels)
           args <- append(list("formula" = defineCase ~ ., "data" = df), args)
           do.call(e1071::naiveBayes, args)
         }, ...)
}

#' Build Linear Discriminant Analysis Model
#'
#' \code{buildLDA} builds a model using the \code{lda} function
#'  from the \code{MASS} package.
#'
#' @inheritParams build.
#' @return Returns an \code{ExprsModel} object.
#' @export
buildLDA <- function(object, top = 0, ...){ # args to lda

  classCheck(object, c("ExprsBinary", "ExprsMulti"),
             "This build method only works for classification tasks.")

  build.(object, top,
         uniqueFx = function(data, labels, ...){

           # Perform linear discriminant analysis via ~ method
           args <- getArgs(...)
           df <- data.frame(data, "defineCase" = labels)
           args <- append(list("formula" = defineCase ~ ., "data" = df), args)
           do.call(MASS::lda, args)
         }, ...)
}

#' Build Support Vector Machine Model
#'
#' \code{buildSVM} builds a model using the \code{svm} function
#'  from the \code{e1071} package.
#'
#' @inheritParams build.
#' @return Returns an \code{ExprsModel} object.
#' @export
buildSVM <- function(object, top = 0, ...){ # args to svm

  classCheck(object, "ExprsArray",
             "This function is applied to the results of ?exprso.")

  build.(object, top,
         uniqueFx = function(data, labels, ...){

           # Perform SVM via ~ method (permits plotting)
           args <- getArgs(...)
           args <- forceArg("probability", TRUE, args)
           args <- forceArg("cross", 0, args)
           df <- data.frame(data, "defineCase" = labels)
           args <- append(list("formula" = defineCase ~ ., "data" = df), args)
           do.call(e1071::svm, args)
         }, ...)
}

#' Build Linear Model
#'
#' \code{buildLM} builds a model using the \code{lm} function.
#'
#' @inheritParams build.
#' @return Returns an \code{ExprsModel} object.
#' @export
buildLM <- function(object, top = 0, ...){ # args to lm

  classCheck(object, "RegrsArray",
             "This feature selection method only works for continuous outcome tasks.")

  build.(object, top,
         uniqueFx = function(data, labels, ...){

           # Perform LM via ~ method
           args <- getArgs(...)
           df <- data.frame(data, "defineCase" = as.numeric(labels))
           args <- append(list("formula" = defineCase ~ ., "data" = df), args)
           do.call(stats::lm, args)
         }, ...)
}

#' Build Generalized Linear Model
#'
#' \code{buildGLM} builds a model using the \code{glm} function.
#'
#' @inheritParams build.
#' @return Returns an \code{ExprsModel} object.
#' @export
buildGLM <- function(object, top = 0, ...){ # args to glm

  classCheck(object, "ExprsArray",
             "This function is applied to the results of ?exprso.")

  if(class(object) == "RegrsArray"){

    build.(object, top,
           uniqueFx = function(data, labels, ...){

             # Perform GLM via ~ method
             args <- getArgs(...)
             df <- data.frame(data, "defineCase" = as.numeric(labels))
             args <- append(list("formula" = defineCase ~ ., "data" = df), args)
             do.call(stats::glm, args)
           }, ...)

  }else if(class(object) %in% c("ExprsBinary", "ExprsMulti")){

    build.(object, top,
           uniqueFx = function(data, labels, ...){

             # Perform GLM via ~ method
             args <- getArgs(...)
             args <- forceArg("family", "binomial", args)
             df <- data.frame(data, "defineCase" = as.numeric(labels) - 1)
             args <- append(list("formula" = defineCase ~ ., "data" = df), args)
             do.call(stats::glm, args)
           }, ...)
  }
}

#' Build Logistic Regression Model
#'
#' \code{buildLR} builds a model using the \code{glm} function.
#'
#' @inheritParams build.
#' @return Returns an \code{ExprsModel} object.
#' @export
buildLR <- function(object, top = 0, ...){ # args to glm

  classCheck(object, c("ExprsBinary", "ExprsMulti"),
             "This build method only works for classification tasks.")

  buildGLM(object, top, ...)
}

#' Build LASSO or Ridge Model
#'
#' \code{buildLASSO} builds a model using the \code{cv.glmnet} function
#'  from the \code{glmnet} package.
#'
#' @inheritParams build.
#' @return Returns an \code{ExprsModel} object.
#' @export
buildLASSO <- function(object, top = 0, ...){ # args to cv.glmnet

  classCheck(object, "ExprsArray",
             "This function is applied to the results of ?exprso.")

  if(class(object) == "RegrsArray"){

    build.(object, top,
           uniqueFx = function(data, labels, ...){

             # Perform LASSO GLM
             args <- getArgs(...)
             args <- defaultArg("nfolds", 5, args)
             args <- append(list("x" = data, "y" = labels), args)
             do.call(glmnet::cv.glmnet, args)
           }, ...)

  }else if(class(object) %in% c("ExprsBinary", "ExprsMulti")){

    build.(object, top,
           uniqueFx = function(data, labels, ...){

             # Perform LASSO GLM
             args <- getArgs(...)
             args <- defaultArg("nfolds", 5, args)
             args <- forceArg("family", "binomial", args)
             args <- append(list("x" = data, "y" = as.numeric(labels) - 1), args)
             do.call(glmnet::cv.glmnet, args)
           }, ...)
  }
}

#' Build Artificial Neural Network Model
#'
#' \code{buildANN} builds a model using the \code{nnet} function
#'  from the \code{nnet} package.
#'
#' @inheritParams build.
#' @return Returns an \code{ExprsModel} object.
#' @export
buildANN <- function(object, top = 0, ...){ # args to nnet

  classCheck(object, "ExprsArray",
             "This function is applied to the results of ?exprso.")

  build.(object, top,
         uniqueFx = function(data, labels, ...){

           # Perform ANN via ~ method
           args <- getArgs(...)
           args <- defaultArg("size", 1, args)
           args <- defaultArg("rang", 1/max(abs(as.vector(data))), args)
           args <- defaultArg("maxit", 1000, args)
           df <- data.frame(data, "defineCase" = labels)
           args <- append(list("formula" = defineCase ~ ., "data" = df), args)
           sink(tempfile())
           m <- do.call(nnet::nnet, args)
           sink()
           m
         }, ...)
}

#' Build Decision Tree Model
#'
#' \code{buildDT} builds a model using the \code{rpart} function
#'  from the \code{rpart} package.
#'
#' Provide \code{cp} as a numeric scalar to trim the \code{rpart} decision tree.
#'  If provided, this argument is passed to the \code{rpart::prune} function.
#'  Set \code{cp = 0} to skip pruning (default behavior).
#'
#' @inheritParams build.
#' @return Returns an \code{ExprsModel} object.
#' @export
buildDT <- function(object, top = 0, ...){ # args to rpart

  classCheck(object, "ExprsArray",
             "This function is applied to the results of ?exprso.")

  build.(object, top,
         uniqueFx = function(data, labels, ...){

           # Perform rpart via ~ method
           args <- getArgs(...)
           args <- defaultArg("cp", 0, args)
           prune <- args$cp
           args <- args[!"cp" == names(args)]

           df <- data.frame(data, "defineCase" = labels)
           args <- append(list("formula" = defineCase ~ ., "data" = df), args)
           tree <- do.call(rpart::rpart, args)
           if(prune != 0) tree <- rpart::prune(tree, cp = prune) # optional pruning
           tree
         }, ...)
}

#' Build Random Forest Model
#'
#' \code{buildRF} builds a model using the \code{randomForest} function
#'  from the \code{randomForest} package.
#'
#' @inheritParams build.
#' @return Returns an \code{ExprsModel} object.
#' @export
buildRF <- function(object, top = 0, ...){ # args to randomForest

  classCheck(object, "ExprsArray",
             "This function is applied to the results of ?exprso.")

  build.(object, top,
         uniqueFx = function(data, labels, ...){

           # Perform RF via ~ method
           args <- getArgs(...)
           df <- data.frame(data, "defineCase" = labels)
           args <- append(list("formula" = defineCase ~ ., "data" = df), args)
           do.call(randomForest::randomForest, args)
         }, ...)
}

#' Build Fuzzy Rule Based Model
#'
#' \code{buildFRB} builds a model using the \code{frbs} function
#'  from the \code{frbs} package.
#'
#' @inheritParams build.
#' @return Returns an \code{ExprsModel} object.
#' @export
buildFRB <- function(object, top = 0, ...){ # args to frbs

  classCheck(object, "ExprsArray",
             "This function is applied to the results of ?exprso.")

  build.(object, top,
         uniqueFx = function(data, labels, ...){

           # Perform frbs via df method
           args <- getArgs(...)
           df <- data.frame(data, "defineCase" = as.numeric(labels))
           args <- append(list("data.train" = df), args)
           do.call(frbs::frbs.learn, args)
         }, ...)
}

#' Build Deep Neural Network Model
#'
#' \code{buildDNN} builds a model using the \code{h2o.deeplearning} function
#'  from the \code{h2o} package.
#'
#' @inheritParams build.
#' @return Returns an \code{ExprsModel} object.
#' @export
buildDNN <- function(object, top = 0, ...){ # args to h2o.deeplearning

  packageCheck("h2o")
  classCheck(object, "ExprsArray",
             "This function is applied to the results of ?exprso.")

  build.(object, top,
         uniqueFx = function(data, labels, ...){

           # Initialize h2o (required)
           localH20 <- h2o::h2o.init()
           args <- getArgs(...)

           # Rename data based on order supplied by top
           colnames(data) <- paste0("id", 1:ncol(data))
           df <- data.frame(data, "defineCase" = as.character(labels))

           # Import data as H2OFrame via a temporary csv
           tempFile <- tempfile(fileext = ".csv")
           utils::write.csv(df, tempFile)
           h2o.data <- h2o::h2o.importFile(path = tempFile, destination_frame = "h2o.data")

           # Prepare arguments and build model
           args <- getArgs(...)
           args <- append(list("x" = colnames(data), "y" = "defineCase",
                               "training_frame" = h2o.data), args)
           do.call(h2o::h2o.deeplearning, args)
         }, ...)
}
