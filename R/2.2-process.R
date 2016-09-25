###########################################################
### Pre-process ExprsArray objects

#' Hard Filter Data
#'
#' \code{modFilter} imposes a hard filter for (gene expression) feature data.
#'
#' This method reproduces the hard filter described by Deb and Reddy (2003)
#'  for pre-processing the hallmark Golub ALL/AML dataset. This filter
#'  first sets all values less than \code{threshold} to \code{threshold}
#'  and all values greater than \code{maximum} to \code{maximum}.
#'
#' Next, this method includes only those features with (a) a range greater
#'  than \code{beta1}, and also (b) a ratio of maximum feature expression to
#'  minimum feature expression greater than \code{beta2}.
#'
#' @param object Specifies the \code{ExprsArray} object to undergo pre-processing.
#' @param threshold A numeric scalar. The value below which to assign this value.
#' @param maximum A numeric scalar. The value above which to assign this value.
#' @param beta1 A numeric scalar. The \code{max - min} range above which to
#'  include the feature. Inclusive with \code{beta2}.
#' @param beta2 A numeric scalar. The \code{max / min} ratio above which to
#'  include the feature. Inclusive with \code{beta1}.
#' @param plotSkip A logical scalar. Toggles whether to produce side-by-side plots
#'  of the before and after \code{ExprsArray} summary.
#'
#' @return A pre-processed \code{ExprsArray} object.
#'
#' @seealso
#' \code{\link{modFilter}}, \code{\link{modTransform}}, \code{\link{modNormalize}}
#'
#' @export
setGeneric("modFilter",
           function(object, threshold, maximum, beta1, beta2, plotSkip = TRUE) standardGeneric("modFilter")
)

#' @describeIn modFilter Method to filter an \code{ExprsArray} object.
#'
#' @export
setMethod("modFilter", "ExprsArray",
          function(object, threshold, maximum, beta1, beta2, plotSkip){

            if(!plotSkip){

              layout(matrix(c(1,4,2,5,3,6), 3, 2, byrow = TRUE))

              # For large datasets, plot only a sub-sample
              v <- as.vector(object@exprs)
              if(length(v) > 1000) v <- sample(v, 1000)
              qqnorm(v, main = "Normal Q-Q Plot (Before)")
              plot(density(v), main = "Density Plot (Before)")
              boxplot(v, horizontal = TRUE, main = "Box Plot (Before)", xlab = "Expression Values")
            }

            # Set minimum and maximum threshold
            object@exprs[object@exprs < threshold] <- threshold
            object@exprs[object@exprs > maximum] <- maximum

            # Calculate RANGE and MAX:MIN
            index.beta1 <- apply(object@exprs, MARGIN = 1, function(gene) max(gene) - min(gene))
            index.beta2 <- apply(object@exprs, MARGIN = 1, function(gene) max(gene) / min(gene))

            # INCLUDE those with variance GREATER THAN beta1 AND beta2
            object@exprs <- object@exprs[index.beta1 > beta1 & index.beta2 > beta2, ]

            if(!plotSkip){

              # For large datasets, plot only a sub-sample
              v <- as.vector(object@exprs)
              if(length(v) > 1000) v <- sample(v, 1000)
              qqnorm(v, main = "Normal Q-Q Plot (After)")
              plot(density(v), main = "Density Plot (After)")
              boxplot(v, horizontal = TRUE, main = "Box Plot (After)", xlab = "Expression Values")

              layout(matrix(c(1), 1, 1, byrow = TRUE))
            }

            return(object)
          }
)

#' Log Transform Data
#'
#' \code{modTransform} log (base 2) transforms feature data.
#'
#' @inheritParams modFilter
#' @return A pre-processed \code{ExprsArray} object.
#'
#' @seealso
#' \code{\link{modFilter}}, \code{\link{modTransform}}, \code{\link{modNormalize}}
#'
#' @export
setGeneric("modTransform",
           function(object, plotSkip = TRUE) standardGeneric("modTransform")
)

#' @describeIn modTransform Method to transform an \code{ExprsArray} object.
#'
#' @export
setMethod("modTransform", "ExprsArray",
          function(object, plotSkip){

            if(!plotSkip){

              layout(matrix(c(1,4,2,5,3,6), 3, 2, byrow = TRUE))

              # For large datasets, plot only a sub-sample
              v <- as.vector(object@exprs)
              if(length(v) > 1000) v <- sample(v, 1000)
              qqnorm(v, main = "Normal Q-Q Plot (Before)")
              plot(density(v), main = "Density Plot (Before)")
              boxplot(v, horizontal = TRUE, main = "Box Plot (Before)", xlab = "Expression Values")
            }

            # Perform log transformation
            object@exprs <- log(object@exprs, base = 2)

            if(!plotSkip){

              # For large datasets, plot only a sub-sample
              v <- as.vector(object@exprs)
              if(length(v) > 1000) v <- sample(v, 1000)
              qqnorm(v, main = "Normal Q-Q Plot (After)")
              plot(density(v), main = "Density Plot (After)")
              boxplot(v, horizontal = TRUE, main = "Box Plot (After)", xlab = "Expression Values")

              layout(matrix(c(1), 1, 1, byrow = TRUE))
            }

            return(object)
          }
)

#' Normalize Data
#'
#' \code{modNormalize} normalizes feature data.
#'
#' This method normalizes subject and/or feature vectors according to the
#'  formula \code{y = (x - mean(x)) / sd(x)}.
#'
#' @inheritParams modFilter
#' @param MARGIN A numeric vector. The margin by which to normalize.
#'  Provide \code{MARGIN = 1} to normalize the feature vector.
#'  Provide \code{MARGIN = 2} to normalize the subject vector.
#'  Provide \code{MARGIN = c(1, 2)} to normalize by the subject vector
#'  and then by the feature vector.
#' @return A pre-processed \code{ExprsArray} object.
#'
#' @seealso
#' \code{\link{modFilter}}, \code{\link{modTransform}}, \code{\link{modNormalize}}
#'
#' @export
setGeneric("modNormalize",
           function(object, MARGIN = c(1, 2), plotSkip = TRUE) standardGeneric("modNormalize")
)

#' @describeIn modNormalize Method to normalize an \code{ExprsArray} object.
#'
#' @export
setMethod("modNormalize", "ExprsArray",
          function(object, MARGIN, plotSkip){

            if(!plotSkip){

              layout(matrix(c(1,4,2,5,3,6), 3, 2, byrow = TRUE))

              # For large datasets, plot only a sub-sample
              v <- as.vector(object@exprs)
              if(length(v) > 1000) v <- sample(v, 1000)
              qqnorm(v, main = "Normal Q-Q Plot (Before)")
              plot(density(v), main = "Density Plot (Before)")
              boxplot(v, horizontal = TRUE, main = "Box Plot (Before)", xlab = "Expression Values")
            }

            # Perform normalization
            if(2 %in% MARGIN) object@exprs <- apply(object@exprs, MARGIN = 2,
                                                    function(samp) (samp - mean(samp)) / sd(samp))

            if(1 %in% MARGIN) object@exprs <- t(apply(object@exprs, MARGIN = 1,
                                                      function(gene) (gene - mean(gene)) / sd(gene)))

            if(!plotSkip){

              # For large datasets, plot only a sub-sample
              v <- as.vector(object@exprs)
              if(length(v) > 1000) v <- sample(v, 1000)
              qqnorm(v, main = "Normal Q-Q Plot (After)")
              plot(density(v), main = "Density Plot (After)")
              boxplot(v, horizontal = TRUE, main = "Box Plot (After)", xlab = "Expression Values")

              layout(matrix(c(1), 1, 1, byrow = TRUE))
            }

            return(object)
          }
)
