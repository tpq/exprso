###########################################################
### Define method for mutating case subjects

#' Swap Case Subjects
#'
#' This experimental function mutates a percentage of case subjects
#'  into noisy positives, false positives, or defined out-groups.
#'
#' This function includes several methods for distorting the features of \code{ExprsBinary}
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

##########################################################
### Cluster subjects by feature data

#' Cluster Subjects
#'
#' This method clusters subjects based on feature data using any one of
#'  seven available clustering algorithms. See Arguments below.
#'
#' Note that this function will expect the argument \code{k} to define the returned
#'  number of clusters, except when \code{how = "kmeans"} in which case this
#'  function will expect the argument \code{centers} instead.
#'
#' @inheritParams fs.
#' @param object An \code{ExprsArray} object. The object containing the subject
#'  data to cluster.
#' @param how A character string. The name of the function used to cluster.
#'  Select from "hclust", "kmeans", "agnes", "clara", "diana", "fanny", or
#'  "pam".
#' @param onlyCluster A logical scalar. Toggles whether to return a processed
#'  cluster object or an updated \code{ExprsArray} object.
#' @param ... Additional arguments to the cluster function and/or
#'  other functions used for clustering (e.g., \code{dist} and
#'  \code{cutree}).
#'
#' @return Typically an \code{ExprsArray} object with subject cluster assignments
#'  added to the \code{$cluster} column of the \code{@@anot} slot.
#'
#' @export
setGeneric("modCluster",
           function(object, top = 0, how = "hclust",
                    onlyCluster = FALSE, ...) standardGeneric("modCluster")
)

#' @describeIn modCluster Method to compare \code{ExprsArray} objects.
#'
#' @importFrom stats hclust kmeans cutree dist
#' @importFrom cluster agnes clara diana fanny pam
#' @export
setMethod("modCluster", "ExprsArray",
          function(object, top, how, onlyCluster, ...){

            args <- as.list(substitute(list(...)))[-1]

            if(class(top) == "numeric"){

              if(length(top) == 1){

                if(top > nrow(object@exprs)) top <- 0
                if(top == 0) top <- nrow(object@exprs)
                top <- rownames(object@exprs[1:top, ])

              }else{

                top <- rownames(object@exprs[top, ])
              }
            }

            data <- t(object@exprs[top, ])

            # Set the argument 'k' aside to use in 'cutree'
            # NOTE: hclust, agnes, and diana use 'k' via 'cutree'
            if(how %in% c("hclust", "agnes", "diana")){

              if(!"k" %in% names(args)){

                cat("Setting 'k' to '2' (default behavior, override explicitly)...\n")
                args <- append(args, list("k" = 2))
              }

              # Store 'k' outside of 'args' list, then remove 'k' argument
              args.k <- args$k
              args <- args[!"k" %in% names(args)]

              if(how == "hclust"){

                # Add 'metric' to args if not already included (i.e., for dist())
                if(!"metric" %in% names(args)){

                  cat("Setting 'metric' to 'euclidean' (default behavior, override explicitly)...\n")
                  args <- append(args, list("metric" = "euclidean"))
                }

                # Add 'p' to args if not already included (i.e., for dist())
                if(!"p" %in% names(args)){

                  cat("Setting 'p' to '2' (default behavior, override explicitly)...\n")
                  args <- append(args, list("p" = 2))
                }

                # Build distance matrix, then remove dist() args
                cat("Calculating dissimilarity matrix using 'metric' and 'p' arguments...\n")
                d <- dist(data, method = args$metric, p = args$p)
                args <- args[!"metric" == names(args)]
                args <- args[!"p" == names(args)]

                args <- append(args, list("d" = d))
                result <- do.call(how, args)

              }else{

                args <- append(args, list("x" = data))
                result <- do.call(how, args)
              }

            }else if(how == "kmeans"){

              # Add 'centers' to args (instead of 'k')
              if(!"centers" %in% names(args)){

                cat("Setting 'centers' to '2' (default behavior, override explicitly)...\n")
                args <- append(args, list("centers" = 2))
              }

              args <- append(args, list("x" = data))
              result <- do.call(how, args)

            }else if(how %in% c("clara", "fanny", "pam")){

              # Add 'k' to args if not already included
              if(!"k" %in% names(args)){

                cat("Setting 'k' to '2' (default behavior, override explicitly)...\n")
                args <- append(args, list("k" = 2))
              }

              args <- append(args, list("x" = data))
              result <- do.call(how, args)

            }else{

              stop("Provided how not recognized. Select from 'hclust', 'kmeans', 'agnes', ",
                   "'clara', 'diana', 'fanny', or 'pam'.\n")
            }

            if(onlyCluster){

              cat("Returning the processed cluster object...\n")
              return(result)

            }else{

              if(how %in% c("hclust", "agnes", "diana")){

                clusts <- cutree(result, k = args.k)
              }

              if(how == "kmeans"){

                clusts <- result$cluster
              }

              if(how %in% c("clara", "fanny", "pam")){

                clusts <- result$clustering
              }

              cat("Updating ExprsArray object...\n")
              object@annot[, "cluster"] <- clusts
              return(object)
            }
          }
)

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
