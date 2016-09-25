###########################################################
### Define method for mutating case subjects

#' Swap Case Subjects
#'
#' This experimental function mutates a percentage of case subjects
#'  into noisy positives, false positives, or defined out-groups.
#'
#' This function includes several methods for distoring the features of \code{ExprsBinary}
#'  subjects. The "rp.1" method randomizes subject vectors to create "subject noise".
#'  The "rp.2" method creates a new subject vector by randomly sampling feature values
#'  from the respective feature vector. The "fp" method creates a new subject vector
#'  by randomly sampling feature values from the respective control feature vector.
#'
#' The "ng" and "tg" methods create out-groups by defining new means for each feature.
#'  These methods yield fixed distributions around new feature means such that
#'  the mean of all new feature means remains constant. The argument \code{theta}
#'  dictates how much the new feature mean might differ from the original feature mean
#'  (where larger \code{theta} values lead to more similar new feature means). For
#'  the "ng" method, the mean of new feature means equals that of the original features
#'  for case subjects only. On the other hand, for the "tg" method, the mean of new
#'  feature means equals that of the original features for all subjects.
#'
#' Alternatively, by providing another \code{ExprsBinary} object as the \code{how}
#'  argument, this function will swap a percentage of case subjects from the main dataset
#'  with control subjects from the second dataset.
#'
#' @param object An \code{ExprsBinary} object to mutate.
#' @param how A character string. The method used to mutate case subjects. Select from
#'  "rp.1", "rp.2", "fp", "ng", or "tg". Alternatively, another \code{ExprsBinary}
#'  object. See Details.
#' @param percent A numeric scalar. The percentage of subjects to mutate.
#' @param theta A numeric scalar. Applies a weight to the distribution of means when
#'  mutating subjects via the "ng" or "tg" method.
#'
#' @return An \code{ExprsBinary} object containing mutated subjects with an index
#'  appended to the \code{$mutated} column of the \code{@@annot} slot.
#'
#' @export
setGeneric("modSwap",
           function(object, how = "fp", percent = 10, theta = 1) standardGeneric("modSwap")
)

#' @describeIn modSwap A method to mutate \code{ExprsBinary} objects.
#'
#' @importFrom stats rnorm
#' @export
setMethod("modSwap", "ExprsBinary",
          function(object, how, percent, theta){

            if(percent < 1 | percent > 100){

              stop("Uh oh! Use an inclusion percentage between 1-100!")
            }

            # Mutate a percent of case subjects
            cases <- object@annot$defineCase %in% "Case"
            mut.size <- round(ncol(object@exprs[, cases]) * percent/100)
            mut.name <- sample(colnames(object@exprs[, cases]), size = mut.size, replace = FALSE)

            # Calculate "before" PCA
            temp1 <- fsPrcomp(object, top = 0)

            if(inherits(how, "ExprsArray")){

              if(sum(how@annot$defineCase %in% "Control") < mut.size){

                stop("Uh oh! Not enough controls in the supplemental 'ExprsBinary' object to fulfill swap!")
              }

              if(!identical(rownames(object@exprs), rownames(how@exprs))){

                stop("Uh oh! The provided 'ExprsBinary' objects do not have matching feature vectors!")
              }

              # Randomly select which controls to use in swap
              from.name <- sample(rownames(how@annot[how@annot$defineCase %in% "Control", ]),
                                  size = mut.size, replace = FALSE)

              # Swap 'object' cases for 'how' controls and store index
              object@exprs[, mut.name] <- how@exprs[, from.name]

            }else if(how == "rp.1"){

              for(mut.col in mut.name){

                object@exprs[, mut.col] <- sample(object@exprs[, mut.col])
              }

            }else if(how == "rp.2"){

              for(mut.col in mut.name){

                object@exprs[, mut.col] <- apply(object@exprs, MARGIN = 1, sample, size = 1)
              }

            }else if(how == "fp"){

              for(mut.col in mut.name){

                object@exprs[, mut.col] <- apply(object@exprs[, !cases], MARGIN = 1, sample, size = 1)
              }

            }else if(how == "ng"){

              # Calculate per-feature case means and sds to make "new group" means and sds
              means <- apply(object@exprs[, cases], MARGIN = 1, mean)
              sds <- apply(object@exprs[, cases], MARGIN = 1, sd)
              ng.means <- unlist(lapply(1:length(means),
                                        function(i) stats::rnorm(1, mean = means[i], sd = sds[i]/theta)))

              for(mut.col in mut.name){

                object@exprs[, mut.col] <-
                  unlist(lapply(1:length(ng.means),
                                function(i) stats::rnorm(1, mean = ng.means[i], sd = sds[i])))
              }

            }else if(how == "tg"){

              # Calculate per-feature means and sds using ALL subjects to make "new group" means and sds
              means <- apply(object@exprs, MARGIN = 1, mean)
              sds <- apply(object@exprs, MARGIN = 1, sd)
              tg.means <- unlist(lapply(1:length(means),
                                        function(i) stats::rnorm(1, mean = means[i], sd = sds[i]/theta)))

              for(mut.col in mut.name){

                object@exprs[, mut.col] <-
                  unlist(lapply(1:length(tg.means),
                                function(i) stats::rnorm(1, mean = tg.means[i], sd = sds[i])))
              }

            }else{

              stop("Provided how not recognized. Select from 'rp.1', 'rp.2', 'fp', 'ng', or 'tg'.")
            }

            # Store Boolean index of mutated subjects in @annot
            object@annot$mutated <- as.numeric(rownames(object@annot) %in% mut.name)

            # Calculate "after" PCA
            temp2 <- fsPrcomp(object, top = 0)

            # Visualize "before" and "after" PCA results
            layout(matrix(c(1,2,3,4), 2, 2, byrow = TRUE))
            plot(temp1, a = 1, b = 2, c = 0, main = "Before", xlab = "PCA1", ylab = "PCA2")
            plot(temp2, a = 1, b = 2, c = 0, main = "After", xlab = "PCA1", ylab = "PCA2")
            plot(temp1, a = 1, b = 3, c = 0, main = "Before", xlab = "PCA1", ylab = "PCA3")
            plot(temp2, a = 1, b = 3, c = 0, main = "After", xlab = "PCA1", ylab = "PCA3")
            layout(matrix(c(1), 1, 1, byrow = TRUE))

            return(object)
          }
)
