##########################################################
### Compare datasets and data subsets

#' Compare \code{ExprsArray} Objects
#'
#' This method compares the values of all \code{ExprsArray} annotations across a
#'  specified annotation term for up to two \code{ExprsArray} objects.
#'  Depending on the composition of each annotation, \code{compare}
#'  will perform either a chi-squared test or an ANOVA test.
#'
#' This method performs two kinds of comparisons. First, it tests all
#'  annotation variables against the annotation supplied by the \code{colBy}
#'  argument for each provided \code{ExprsArray} object. In other words,
#'  the \code{colBy} argument determines which annotation to use as the
#'  independent variable for "internal" comparisons. Second, it tests
#'  all annotation variables between the provided \code{ExprsArray} objects.
#'  Providing \code{array.valid = NULL} will skip the between comparisons.
#'
#' This method will test annotations using either a chi-squared test or an
#'  ANOVA test depending on the class of the values stored by the tested column.
#'  The presence of a "character" or "factor" in the tested column will trigger
#'  a chi-squared test. As such, this method requires the user to select
#'  a \code{colBy} annotation that contains categorical data (i.e., to use as
#'  the independent variable).
#'
#' We anticipate that this method will serve as a useful adjunct to
#'  \code{\link{modCluster}}. However, it may also help in quickly determining
#'  whether the data \code{\link{split}} has yielded comparable training and
#'  test sets in terms of the annotations included in \code{@@annot}.
#'
#' @param object The \code{ExprsArray} object used when comparing annotations.
#' @param array.valid A second \code{ExprsArray} object used when comparing
#'  annotations. Optional. Exclude with \code{array.valid = NULL}.
#' @param colBy A character string. The annotation column against which to
#'  compare all other annotation terms (i.e., to test as the independent
#'  variable).
#' @param cutoff A numeric scalar. The p-value cutoff that determines when
#'  the annotation test returns a \code{TRUE} result
#'
#' @return A list of three logical vectors. The first and second elements
#'  of the list correspond to "internal" comparisons for the two provided
#'  \code{ExprsArray} objects, respectively. The third element of the list
#'  corresponds to comparisons made between the provided objects.
#'
#' @export
setGeneric("compare",
           function(object, array.valid = NULL, colBy = "defineCase",
                    cutoff = .05) standardGeneric("compare")
)

#' @describeIn compare Method to compare \code{ExprsArray} objects.
#'
#' @importFrom stats na.omit lm anova chisq.test
#' @export
setMethod("compare", "ExprsArray",
          function(object, array.valid, colBy, cutoff){

            if(length(colBy) > 1){

              stop("Uh oh! You can only include one 'colBy' annotation term.")
            }

            if(!is.null(array.valid)){

              if(!inherits(object, "ExprsArray")){

                stop("Uh oh! You can only compare ExprsArray objects.")
              }
            }

            # Perform ANOVA or chi-square across each provided array
            results <- vector("list", 3)
            for(i in 1:length(results)){

              if(i == 1){

                cat("## Making internal comparisons for first object (across 'colBy' column)...\n", sep = "")
                test.data <- object@annot
              }

              if(i == 2 & !is.null(array.valid)){

                cat("## Making internal comparisons for second object (across 'colBy' column)...\n", sep = "")
                test.data <- array.valid@annot
              }

              if(i == 3 & !is.null(array.valid)){

                cat("## Making comparisons between objects (independent of 'colBy' column)...\n", sep = "")
                test.data <- rbind(data.frame("array" = "object",
                                              object@annot),
                                   data.frame("array" = "array.valid",
                                              array.valid@annot[, colnames(object@annot)]))

                # For the third comparison, make 'colBy' represent the data source
                # This forces a comparison between 'object' and 'array.valid'
                colBy <- "array"
              }

              # Make sure colBy contains categorical data for both 'object' and 'array.valid'
              if(!class(test.data[, colBy]) %in% c("character", "factor")){

                stop("Uh oh! The 'colBy' variable must contain categorical data.")
              }

              # Check whether the values of 'annot' differ across 'colBy' for each data source
              # Then, check whether the values of 'annot' differ across the data sources
              if(i == 1 | !is.null(array.valid)){

                # Prepares results vector only if running checks
                annots <- colnames(test.data)[!colnames(test.data) %in% colBy]
                results[[i]] <- vector("logical", length(annots))
                names(results[[i]]) <- annots

                for(annot in annots){

                  if(class(test.data[, annot]) %in% c("character", "factor")){

                    # Perform Chi-squared test and save significance as boolean
                    compare <- table(test.data[, c(annot, colBy)])
                    chisq.p <- chisq.test(compare)$p.value
                    cat("\t", "Annotation ", annot, ": ", chisq.p, "\n", sep = "")
                    results[[i]][annot] <- ifelse(!is.na(chisq.p), chisq.p < cutoff, FALSE)

                  }else{

                    df <- data.frame("y" = test.data[, annot],
                                     "x" = as.factor(test.data[, colBy]))

                    # Check if ANOVA will work after removing NA values
                    if(length(unique(stats::na.omit(df)$x)) > 1){

                      # Fit linear model to independent variable
                      fit <- lm(y ~ x, data = df, na.action = stats::na.omit)

                      # Perform ANOVA test and save significance as boolean
                      anova.p <- anova(fit)$Pr[1]
                      cat("\t", "Annotation ", annot, ": ", anova.p, "\n", sep = "")
                      results[[i]][annot] <- anova.p < cutoff

                    }else{

                      cat("\t", "Annotation ", annot, ": insufficient measurements\n", sep = "")
                      results[[i]][annot] <- FALSE
                    }
                  }
                }
              }
            }

            return(results)
          }
)
