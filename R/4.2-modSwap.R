###########################################################
### Define method for swapping cases with controls

#' Swap Case Subjects
#'
#' This experimental function swaps a percentage of case subjects from one dataset with
#'  control subjects from another dataset.
#'
#' @param object An \code{ExprsBinary} object containing cases to swap with controls.
#'
#' @export
setGeneric("modSwap",
           function(object, ...) standardGeneric("modSwap")
)

#' @describeIn modSwap A method to swap \code{ExprsBinary} objects.
#'
#' @param from An \code{ExprsBinary} object from which to draw the control subjects.
#' @param percent A numeric scalar. The percentage of subjects to mutate.
#'
#' @return An \code{ExprsBinary} object containing swapped subjects with an index
#'  appended to the \code{$mutated} column of the \code{@@annot} slot.
#'
#' @export
setMethod("modSwap", "ExprsBinary",
          function(object, from, percent = 10){

            if(percent < 1 | percent > 100){

              stop("You must choose an inclusion percentage between 1-100!")
            }

            if(class(from) != "ExprsBinary"){

              stop("You must provide an 'ExprsBinary' object for the 'from' argument!")
            }

            # Swap a percent of case subjects
            cases <- object@annot$defineCase %in% "Case"
            mut.size <- round(ncol(object@exprs[, cases]) * percent/100)
            mut.name <- sample(colnames(object@exprs[, cases]), size = mut.size, replace = FALSE)

            # Calculate "before" PCA
            temp1 <- fsPrcomp(object, probes = 0)

            if(sum(from@annot$defineCase %in% "Control") < mut.size){

              stop("Uh oh! Not enough controls in the supplemental 'ExprsBinary' object to fulfill swap!")
            }

            if(!identical(rownames(object@exprs), rownames(from@exprs))){

              stop("Uh oh! The provided 'ExprsBinary' objects do not have matching probe vectors!")
            }

            # Randomly select which controls to use in swap
            from.name <- sample(rownames(from@annot[from@annot$defineCase %in% "Control", ]),
                                size = mut.size, replace = FALSE)

            # Swap 'object' cases for 'from' controls and store index
            object@exprs[, mut.name] <- from@exprs[, from.name]
            object@annot$mutated <- as.factor(rownames(object@annot) %in% mut.name)

            # Calculate "after" PCA
            temp2 <- fsPrcomp(object, probes = 0)

            # Visualize "before" and "after" PCA results
            layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
            plot(temp1@exprs[1, ], temp1@exprs[2, ], col = ifelse(cases, "red", "black"),
                 pch = ifelse(object@annot$mutated, 22, 16), main = "Before", xlab = "PCA1", ylab = "PCA2")
            plot(temp2@exprs[1, ], temp2@exprs[2, ], col = ifelse(cases, "red", "black"),
                 pch = ifelse(object@annot$mutated, 22, 16), main = "After", xlab = "PCA1", ylab = "PCA2")
            plot(temp1@exprs[1, ], temp1@exprs[3, ], col = ifelse(cases, "red", "black"),
                 pch = ifelse(object@annot$mutated, 22, 16), main = "Before", xlab = "PCA1", ylab = "PCA3")
            plot(temp2@exprs[1, ], temp2@exprs[3, ], col = ifelse(cases, "red", "black"),
                 pch = ifelse(object@annot$mutated, 22, 16), main = "After", xlab = "PCA1", ylab = "PCA3")
            layout(matrix(c(1), 1, 1, byrow = TRUE))

            return(object)
          }
)
