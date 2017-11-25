###########################################################
### Import data from data.frame, eSet, or file

#' Import Data as ExprsArray
#'
#' A convenience function that builds an \code{ExprsArray} object.
#'  This function is no longer supported. Please use \code{\link{exprso}} instead.
#'
#' Importing a \code{data.frame} object:
#'
#' This function expects that the imported \code{data.frame} has the following format:
#'  rows indicate subject entries while columns indicate measured variables.
#'  The first several columns should contain annotation information (e.g., age, sex, diagnosis).
#'  The remaining columns should contain feature data (e.g., expression values).
#'  The argument \code{begin} defines the j-th column at which the feature
#'  data starts. This function automatically removes any features with \code{NA} values.
#'  Take care to remove any \code{factor} columns before importing.
#'
#' Importing an \code{ExpressionSet} object:
#'
#' The package Biobase maintains a popular class object called \code{ExpressionSet} that
#'  often gets used to store expression data. This function converts this \code{eSet}
#'  object into an \code{ExprsArray} object. This function automatically removes any
#'  features with \code{NA} values.
#'
#' Importing a \code{file}:
#'
#' \code{arrayExprs} can also build an \code{ExprsArray} object from a tab-delimited
#'  data file, passing along the \code{file} and \code{...} argument(s) to
#'  \code{\link{read.delim}}. All rules for \code{data.frame} import also apply here.
#'  By default, \code{arrayExprs} forces \code{stringsAsFactors = FASE}.
#'
#' @param object What to import as an \code{ExprsArray} object. See Details.
#' @param colBy A numeric or character index. The column that contains group annotations.
#' @param include A list of character vectors. Specifies which annotations in \code{colBy}
#'  to include in which groups. Each element of the list specifies a unique group while
#'  each element of the character vector specifies which annotations define that group. For
#'  binary classification, the first list element defines the negative, or control, group.
#' @param colID A numeric or character index. The column used to name subjects.
#'  For \code{data.frame} or file import only.
#' @param begin A numeric scalar. The j-th column at which feature data starts.
#'  For \code{data.frame} or file import only.
#' @param ... Additional arguments passed along to \code{read.delim}.
#'  For file import only.
#'
#' @return An \code{ExprsArray} object.
#'
#' @seealso
#' \code{\link{ExprsArray-class}}, \code{\link{GSE2eSet}}
#' @importFrom utils read.delim
#' @export
arrayExprs <- function(object, colBy, include, colID, begin, ...){

  if(missing(colBy)) stop("Uh oh! User must specify 'colBy' argument!")
  if(missing(include)) stop("Uh oh! User must specify 'include' argument!")
  if(!class(include) == "list") stop("Uh oh! User must provide 'include' argument as list!")
  if(length(include) < 2) stop("Uh oh! User must provide at least two classes!")

  if(class(object) == "data.frame"){

    if(missing(colID)) stop("Uh oh! User must specify 'colID' argument!")
    if(missing(begin)) stop("Uh oh! User must specify 'begin' argument!")
    rownames(object) <- object[, colID]
    exprs <- exprs <- t(object[, begin:ncol(object)])
    annot <- object[, 1:(begin-1)]

  }else if(class(object) == "ExpressionSet"){

    packageCheck("Biobase")

    exprs <- Biobase::exprs(object)
    annot <- object@phenoData@data

  }else if(class(object) == "character" & file.exists(object)){

    if(missing(colID)) stop("Uh oh! User must specify 'colID' argument!")
    if(missing(begin)) stop("Uh oh! User must specify 'begin' argument!")
    args <- getArgs(...)
    args <- forceArg("stringsAsFactors", FALSE, args)
    args <- append(args, list("file" = object))
    object <- do.call("read.delim", args)
    rownames(object) <- object[, colID]
    exprs <- t(object[, begin:ncol(object)])
    annot <- object[, 1:(begin-1)]

  }else{

    stop("Uh oh! No default method for importing this object as an ExprsArray.")
  }

  newClass <- ifelse(length(include) == 2, "ExprsBinary", "ExprsMulti")
  array <- new(newClass, exprs = exprs, annot = annot, preFilter = NULL, reductionModel = NULL)

  # Force @annot rownames to mirror proper @exprs colnames
  colnames(array@exprs) <- make.names(colnames(array@exprs), unique = TRUE)
  rownames(array@annot) <- colnames(array@exprs)

  # Use 'include' to filter subjects and label classes
  array@annot <- array@annot[array@annot[, colBy] %in% unlist(include), ]
  array@exprs <- array@exprs[, rownames(array@annot)]
  for(i in 1:length(include)){

    array@annot[array@annot[, colBy] %in% include[[i]], "defineCase"] <- i
  }

  if(newClass == "ExprsBinary"){

    array@annot$defineCase <- ifelse(array@annot$defineCase == 1, "Control", "Case")

  }else if(newClass == "ExprsMulti"){

    array@annot$defineCase <- factor(array@annot$defineCase)

  } # if regrso, do nothing here

  # Remove features with missing values
  if(any(is.na(array@exprs))){

    cat("Removing features with missing values...\n")
    array@exprs <- array@exprs[apply(array@exprs, 1, function(x) !any(is.na(x))), ]
  }

  return(array)
}

#' Convert GSE to eSet
#'
#' A convenience function that builds an \code{eSet} object from a GSE data source.
#'
#' The NCBI GEO hosts files in GSE or GDS format, the latter of which exists as a curated version
#'  the former. These GDS data files easily convert to an  \code{ExpressionSet} (abbreviated
#'  \code{eSet}) object using the \code{GDS2eSet} function available from the GEOquery package.
#'  However, not all GSE data files have a corresponding GDS data file available. To convert GSE
#'  data files into \code{eSet} objects, \code{exprso} provides this convenience function.
#'
#' However, the user should note that GSE data files do not always get stored in an easy to parse format.
#'  Although this function has worked successfully with some GSE data files, we cannot make any
#'  guarantee that it will work for all GSE data files.
#'
#' To acquire GSE data files, use the function \code{getGEO} from the GEOquery package (e.g.,
#'  \code{getGEO("GSExxxxx", GSEMatrix = FALSE)}). For more information, see the GEOquery package.
#'
#' @param gse A GSE data object retrieved using GEOquery.
#' @param colBy A character string. The GSE column name that contains the feature value.
#'  If missing, function will prompt user for a column name after previewing options.
#' @param colID A character string. The GSE column name that contains the feature identity.
#'  If missing, function will prompt user for a column name after previewing options.
#'
#' @return An \code{ExpressionSet} object.
#'
#' @seealso
#' \code{\link{ExprsArray-class}}, \code{\link{arrayExprs}}
#' @export
GSE2eSet <- function(gse, colBy, colID){

  packageCheck("GEOquery")
  packageCheck("Biobase")
  packageCheck("affy")

  # Check for non-unique platforms
  gsms <- unlist(lapply(GEOquery::GSMList(gse), function(g){ GEOquery::Meta(g)$platform}))
  if(length(unique(gsms)) > 1) stop("GSE contains non-unique platforms!")

  # Provide an opportunity for user to select a new platform ID column
  if(missing(colID)){

    cat("The columns available for platform ID include:\n")
    print(GEOquery::Columns(GEOquery::GSMList(gse)[[1]]))
    cat("\n")
    colID <- readline(prompt = "Which column (by name) will you use for platform ID? \n")
  }

  # Provide an opportunity for user to select a new platform VALUE column
  if(missing(colBy)){

    cat("The columns available for platform VALUE include:\n")
    print(GEOquery::Columns(GEOquery::GSMList(gse)[[1]]))
    cat("\n")
    colBy <- readline(prompt = "Which column (by name) will you use for platform VALUE? \n")
  }

  # Establish a template for all features
  featsets <- GEOquery::Table(GEOquery::GPLList(gse)[[1]])$ID

  # Retrieve feature values for each sample
  vals <- lapply(GEOquery::GSMList(gse),
                 function(g){

                   as.numeric(GEOquery::Table(g)[match(featsets, GEOquery::Table(g)[, colID]), colBy])
                 }
  )

  # Combine samples into data.frame
  exprs <- data.frame(vals, row.names = featsets)

  # Retrieve annotations for each sample
  pdata <- lapply(GEOquery::GSMList(gse),
                  function(g){

                    # Retrieve all columns containing "characteristics_ch" in name
                    annotIndex <- grepl("characteristics_ch", names(GEOquery::Meta(g)))
                    annotSlots <- names(GEOquery::Meta(g))[annotIndex]

                    # If there is comma delimitation and multiple annotations
                    if(any(grepl(",", unlist(GEOquery::Meta(g)[annotSlots])) &
                           unlist(lapply(GEOquery::Meta(g)[annotSlots], function(x) length(x) == 1)))){

                      # Split comma containing annotations
                      characteristics <- unlist(strsplit(unlist(GEOquery::Meta(g)[annotSlots]), ",\\s*"))

                    }else{

                      # Pass along annotations unsplit
                      characteristics <- unlist(GEOquery::Meta(g)[annotSlots])
                    }

                    # Retrieve substance of "characteristics" without characteristic name
                    annots <- lapply(strsplit(characteristics, ":\\s*"), function(split) split[2])

                    # Prepare annotations in AnnotatedDataFrame format
                    annots <- data.frame(annots, row.names = GEOquery::Meta(g)$geo_accession)

                    # Rename columns based on "characteristics" name
                    colnames(annots) <- unlist(lapply(strsplit(characteristics, ":\\s*"), function(x) x[1]))

                    # Add a "sample" column
                    cbind(data.frame("sample" = GEOquery::Meta(g)$geo_accession), annots)
                  }
  )

  # Combine annotations into data.frame
  phenoData <- do.call(plyr::rbind.fill, pdata)
  rownames(phenoData) <- phenoData$sample

  # Build eSet object
  eset <- new("ExpressionSet", exprs = as.matrix(exprs), phenoData = as(phenoData, "AnnotatedDataFrame"))

  return(eset)
}

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

##########################################################
### Generalize multi-class feature selection

#' Serialize "1 vs. all" Feature Selection
#'
#' This experimental function converts multiple feature rank lists,
#'  derived from "1 vs. all" binary feature selection, into a single
#'  feature rank list. This function is not in use in this package.
#'
#' After passing a feature selection method through \code{doMulti},
#'  a set of ranked features gets returned for each one of the
#'  total number of levels in the \code{$defineCase} factor. In
#'  order to proceed with model deployment (at least in the setting
#'  of a conventional pipeline where feature selection occurs
#'  prior to classifier construction), multiple feature rankings
#'  would need to get serialized into a single feature rank list.
#'  \code{reRank} accomplishes this by calculating the rank sum
#'  for each feature across all "1 vs. all" feature selection
#'  tasks. Features found in one rank list but not in another
#'  receive a numeric rank equal to one more than the maximum rank
#'  in that feature rank list. The presence of a NA placeholder
#'  (see: \code{\link{doMulti}}) will not impact \code{reRank}.
#'
#' We note here, however, that a better approach would deploy
#'  "1 vs. all" feature selection and classifier construction
#'  simultaneously, rather than "1 vs. all" feature selection
#'  followed by "1 vs. all" classifier construction. This is
#'  now implemented as \code{\link{plGridMulti}}.
#'
#' @param fss The result of a \code{doMulti} function call.
#' @return A vector of re-ranked features. See Details.
#' @export
reRank <- function(fss){

  # Remove any NULL placeholders
  fss <- fss[!sapply(fss, is.null)]

  # Line up rankings for each doMulti result
  i <- 1
  while(i < length(fss)){

    if(i == 1){

      rank <- data.frame("feat" = rownames(fss[[i]]@exprs),
                         "rank" = 1:nrow(fss[[i]]@exprs))
    }

    rank.next <- data.frame("feat" = rownames(fss[[i + 1]]@exprs),
                            "rank" = 1:nrow(fss[[i + 1]]@exprs))

    rank <- merge(rank, rank.next,
                  by.x = "feat",
                  by.y = "feat",
                  all = TRUE)

    i <- i + 1
  }

  # Clean up rank table
  rownames(rank) <- rank[, "feat"]
  rank <- rank[, !colnames(rank) %in% "feat"]

  # For each class, replace any NAs with 1 more than the maximum rank
  for(col in 1:ncol(rank)){

    maximum <- max(rank[!is.na(rank[, col]), col])
    rank[is.na(rank[, col]), col] <- maximum + 1
  }

  # Add per-class ranks to make final rank list
  final <- rowSums(rank)
  final <- names(final)[order(final)]

  return(final)
}
