###########################################################
### Define generic functions

setGeneric("arraySubset",
           function(object, ...) standardGeneric("arraySubset")
)

###########################################################
### Subset data

setMethod("arraySubset", "ExprsArray",
          function(object, colBy, include){

            # Filter out samples not matching set.include
            object@annot <- object@annot[object@annot[, colBy] %in% include, ]

            # Subset @exprs to filter samples no longer provided by @annot
            object@exprs <- object@exprs[, rownames(object@annot)]

            return(object)
          }
)
