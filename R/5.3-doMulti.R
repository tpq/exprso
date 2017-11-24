#' Perform "1 vs. all" Task
#'
#' A function to execute multiple "1 vs. all" binary tasks.
#'
#' \code{doMulti} depends on the total number of levels in the
#'  \code{$defineCase} factor. If a training set is missing any
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
#' @inheritParams fs.
#' @param method A character string. The method to apply.
#' @return A list of the results from \code{method}.
#' @export
setGeneric("doMulti",
           function(object, top, method, ...) standardGeneric("doMulti")
)

#' @describeIn doMulti Method to execute multiple "1 vs. all" binary tasks.
#' @export
setMethod("doMulti", "ExprsMulti",
          function(object, top, method, ...){

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
)
