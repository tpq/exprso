###########################################################
### Define generic functions

setGeneric("splitSample",
           function(object, ...) standardGeneric("splitSample")
)

setGeneric("splitStratify",
           function(object, ...) standardGeneric("splitStratify")
)

###########################################################
### Split data

# Set replace = FALSE to sample training set WITHOUT replacement
# Set percent.include = 100 to return a NULL demi-holdout array
setMethod("splitSample", "ExprsArray",
          function(object, percent.include, ...){ # args to sample

            warning("This method is not truly random; at least one of every class will appear in validation set!")

            if(percent.include < 1 | percent.include > 100) stop("Uh oh! Use an inclusion percentage between 1-100!")

            # Calculate size of bootstrap set
            size <- round((ncol(object@exprs) * percent.include)/100, digits = 0)

            all.in <- FALSE # Set all.in to FALSE by default
            counter <- 1
            while(!all.in){ # Until sample has at least one of every class

              # Terminate after 10 iterations
              counter <- counter + 1
              if(counter > 10) stop("splitSample could not find a solution. Check the supplied parameters.")

              # Perform a simple random sample with or without replacement
              boot <- sample(colnames(object@exprs), size = size, ...)

              # Check if sample has at least one of every class
              if(all(unique(object@annot$defineCase) %in% unique(object@annot[boot, "defineCase"]))) all.in <- TRUE
            }

            # Build bootstrap set
            cat("Building the random training set...\n\n")
            array.train <- new(class(object),
                               exprs = as.matrix(object@exprs[, boot]),
                               annot = object@annot[boot, ],
                               preFilter = object@preFilter,
                               reductionModel = object@reductionModel)

            # Clean up 1-subject artifact
            colnames(array.train@exprs) <- rownames(array.train@annot)

            # Build demi-holdout set based on the size of the bootstrap set
            if(ncol(array.train@exprs) < ncol(object@exprs)){

              cat("Building a validation set with remaining samples...\n\n")
              array.valid <- new(class(object),
                                 exprs = as.matrix(object@exprs[, !colnames(object@exprs) %in% boot]),
                                 annot = object@annot[!rownames(object@annot) %in% boot, ],
                                 preFilter = object@preFilter,
                                 reductionModel = object@reductionModel)

              # Clean up 1-subject artifact
              colnames(array.valid@exprs) <- rownames(array.valid@annot)

            }else{

              cat("Building a NULL validation set...\n\n")
              array.valid <- NULL
            }

            return(list("array.train" = array.train, "array.valid" = array.valid))
          }
)

#'
#'
#'
#'
#'
#'
#' @import sampling
#' @export
setMethod("splitStratify", "ExprsBinary",
          function(object, percent.include, colBy = NULL, bin, breaks, ...){ # args to cut

            require(sampling)

            if(percent.include < 1 | percent.include > 100) stop("Uh oh! Use an inclusion percentage between 1-100!")

            if(!is.null(colBy)){

              # Perform pre-stratification binning
              vals <- lapply(seq_along(bin),
                             function(i){

                               # For each colBy argument, bin when appropriate
                               if(bin[i]) return(cut(object@annot[, colBy[i]], breaks = breaks[[i]], ...))
                               if(!bin[i]) return(object@annot[, colBy[i]])
                             }
              )

              # Add "defineCase" values to the vals bin list
              vals <- c(list(object@annot[, "defineCase"]), vals)

              # Name vals according to source
              names(vals) <- c("defineCase", colBy)

              # Build a sampling data.frame
              df <- data.frame(vals,
                               row.names = rownames(object@annot),
                               stringsAsFactors = FALSE)

              # Order vals
              index <- do.call(order, vals)
              df <- df[index, ]

              # Remove any row with NAs
              index.NAs <- apply(df, 1, function(row) !NA %in% row)
              df <- df[index.NAs, ]

              cat("\nPre-stratification table:\n")
              print(table(df))

              # Manipulate order of apply(df) so that size vector matches strata expectations
              sizes <- apply(table(df[, c("defineCase", rev(colBy))]), MARGIN = -1, FUN = min)
              sizes <- round(sizes * percent.include/100)
              sizes <- rep(sizes, 2)

              # Provide error for anticipated stratum of size 0
              if(0 %in% sizes) stop("This function cannot create a stratum of size 0. Subset first!")

              # Perform stratification with remaining non-NA colBy values (e.g. those introduced by binning)
              s <- sampling::strata(df, stratanames = colnames(df), size = sizes, method = "srswor")
              if(!identical(rownames(s), rownames(df)[s$ID_unit])) stop("Uh-oh! DEBUG ERROR: 001")

              cat("\nWeighted stratification results:\n")
              print(table(s[, colnames(df)]))
            }

            if(is.null(colBy)){

              # Build a sampling data.frame
              df <- data.frame("defineCase" = object@annot$defineCase,
                               row.names = rownames(object@annot),
                               stringsAsFactors = FALSE)

              cat("\nPre-stratification table:\n")
              print(table(df))

              # Compute strata sizes
              sizes <- min(table(df))
              sizes <- round(sizes * percent.include/100)
              sizes <- rep(sizes, 2)

              # Provide error for anticipated stratum of size 0
              if(0 %in% sizes) stop("This function cannot create a stratum of size 0. Subset first!")

              # Stratify
              s <- sampling::strata(df, stratanames = colnames(df), size = sizes, method = "srswor")
              rownames(s) <- rownames(df)[s$ID_unit]

              cat("\nWeighted stratification results:\n")
              print(table(s[, colnames(df)]))
            }

            # Build array.train from the stratified sample
            cat("\nBuilding the stratified training set...\n\n")
            array.train <- new(class(object),
                               exprs = object@exprs[, rownames(s)],
                               annot = object@annot[rownames(s), ],
                               preFilter = object@preFilter,
                               reductionModel = object@reductionModel)

            # Dump remaining samples into array.valid
            cat("Building a validation set with remaining samples...\n\n")
            array.valid <- new(class(object),
                               exprs = object@exprs[, !colnames(object@exprs) %in% colnames(array.train@exprs)],
                               annot = object@annot[!rownames(object@annot) %in% rownames(array.train@annot), ],
                               preFilter = object@preFilter,
                               reductionModel = object@reductionModel)

            return(list("array.train" = array.train, "array.valid" = array.valid))
          }
)
