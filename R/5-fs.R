###########################################################
### Define generic functions

setGeneric("fsSample", function(object, ...) standardGeneric("fsSample"))
setGeneric("fsStats", function(object, ...) standardGeneric("fsStats"))
setGeneric("fsPrcomp", function(object, ...) standardGeneric("fsPrcomp"))
setGeneric("fsPenalizedSVM", function(object, ...) standardGeneric("fsPenalizedSVM"))
setGeneric("fsPathClassRFE", function(object, ...) standardGeneric("fsPathClassRFE"))
setGeneric("fsEbayes", function(object, ...) standardGeneric("fsEbayes"))
setGeneric("fsMrmre", function(object, ...) standardGeneric("fsMrmre"))

###########################################################
### Select features

# NOTE: @preFilter stores feature selection history and exists to recreate fs during svmPredict
# NOTE: IF probes = 0, include ALL probes when building data; otherwise, select top N occurring first
# NOTE: IF probes is a character vector, include only these probes when building data
# NOTE: 'probes' argument refers to what you feed INTO the function, not what you expect OUT

setMethod("fsSample", "ExprsBinary",
          function(object, probes, ...){ #args to ebayes

            # Convert 'numeric' probe argument to 'character' probe vector
            if(class(probes) == "numeric"){

              if(probes == 0) probes <- nrow(object@exprs)
              probes <- rownames(object@exprs[1:probes, ])
            }

            # Build data using supplied 'character' probe vector
            if(class(probes) == "character"){

              data <- t(object@exprs[probes, ])
            }

            # Randomly sample probes
            final <- sample(probes)

            array <- new("ExprsBinary",
                         exprs = object@exprs[final,],
                         annot = object@annot,
                         preFilter = append(object@preFilter, list(final)),
                         reductionModel = append(object@reductionModel, list(NA))
            )

            return(array)
          }
)

setMethod("fsStats", "ExprsBinary",
          function(object, probes, how = "t.test", ...){ # args to ks.test, ks.boot, or t.test

            # Convert 'numeric' probe argument to 'character' probe vector
            if(class(probes) == "numeric"){

              if(probes == 0) probes <- nrow(object@exprs)
              probes <- rownames(object@exprs[1:probes, ])
              data <- t(object@exprs[probes, ])
            }

            # Build data using supplied 'character' probe vector
            if(class(probes) == "character"){

              data <- t(object@exprs[probes, ])
            }

            # Retrieve case and control subjectIDs
            cases <- object@annot$defineCase %in% "Case"
            conts <- object@annot$defineCase %in% "Control"

            # Load dependencies for ks.boot
            if(how == "ks.boot") require(Matching)

            # Initialize p-value container
            p <- vector("numeric", length(probes))

            for(i in 1:length(probes)){

              if(how == "ks.test"){ p[i] <- ks.test(object@exprs[i, cases], object@exprs[i, conts], ...)$p.value
              }else if(how == "ks.boot"){ p[i] <- ks.boot(object@exprs[i, cases], object@exprs[i, conts], ...)$ks.boot.pvalue
              }else if(how == "t.test"){ p[i] <- t.test(object@exprs[i, cases], object@exprs[i, conts], ...)$p.value
              }else{ stop("Uh oh! Provided 'how' argument not recognized!")}
            }

            final <- probes[order(p)]
            array <- new("ExprsBinary",
                         exprs = object@exprs[final,],
                         annot = object@annot,
                         preFilter = append(object@preFilter, list(final)),
                         reductionModel = append(object@reductionModel, list(NA))
            )

            return(array)
          }
)

setMethod("fsPrcomp", "ExprsBinary",
          function(object, probes, ...){ # args to prcomp

            # Convert 'numeric' probe argument to 'character' probe vector
            if(class(probes) == "numeric"){

              if(probes == 0) probes <- nrow(object@exprs)
              probes <- rownames(object@exprs[1:probes, ])
              data <- t(object@exprs[probes, ])
            }

            # Build data using supplied 'character' probe vector
            if(class(probes) == "character"){

              data <- t(object@exprs[probes, ])
            }

            # NOTE: as.data.frame will not rename columns
            data <- data.frame(data)

            # ATTENTION: We want dependent variables as columns
            reductionModel <- prcomp(data, ...)

            cat("\nDimension reduction model summary:\n\n")
            print(summary(reductionModel))

            # ATTENTION: The value of predict(reductionModel, data) equals $x
            # @preFilter stores probes used to build reductionModel (i.e. as passed on by 'probes' argument)
            # This information will automatically distill the data when calling modHistory
            array <- new("ExprsBinary",
                         exprs = t(reductionModel$x),
                         annot = object@annot,
                         preFilter = append(object@preFilter, list(probes)),
                         reductionModel = append(object@reductionModel, list(reductionModel))
            )

            return(array)
          }
)

setMethod("fsPenalizedSVM", "ExprsBinary",
          function(object, probes, ...){ # args to svm.fs

            require(penalizedSVM)

            # Convert 'numeric' probe argument to 'character' probe vector
            if(class(probes) == "numeric"){

              if(probes == 0) probes <- nrow(object@exprs)
              probes <- rownames(object@exprs[1:probes, ])
              data <- t(object@exprs[probes, ])
            }

            # Build data using supplied 'character' probe vector
            if(class(probes) == "character"){

              data <- t(object@exprs[probes, ])
            }

            # Set labels as factor
            labels <- factor(object@annot[rownames(data), "defineCase"], levels = c("Control", "Case"))
            levels(labels) <- c(-1, 1)

            # Run svm.fs()
            fs <- svm.fs(x = data, y = labels, ...)

            final <- names(fs$model$w)

            if(length(final) < 2) stop("Uh oh! fsPenalizedSVM did not find enough features!")
            array <- new("ExprsBinary",
                         exprs = object@exprs[final, ],
                         annot = object@annot,
                         preFilter = append(object@preFilter, list(final)),
                         reductionModel = append(object@reductionModel, list(NA))
            )

            return(array)
          }
)

setMethod("fsPathClassRFE", "ExprsBinary",
          function(object, probes, ...){ # args to fit.rfe

            require(pathClass)

            # Convert 'numeric' probe argument to 'character' probe vector
            if(class(probes) == "numeric"){

              if(probes == 0) probes <- nrow(object@exprs)
              probes <- rownames(object@exprs[1:probes, ])
              data <- t(object@exprs[probes, ])
            }

            # Build data using supplied 'character' probe vector
            if(class(probes) == "character"){

              data <- t(object@exprs[probes, ])
            }

            # Set labels as factor
            labels <- factor(object@annot[rownames(data), "defineCase"], levels = c("Control", "Case"))

            # Set up "make.names" key for improper @exprs row.names
            key <- data.frame("old" = colnames(data), "new" = make.names(colnames(data)))

            # NOTE: RFE as assembled by pathClass is via a linear kernel only
            # NOTE: By default, fit.rfe iterates through C = 10^c(-3:3)
            # Run fit.rfe()
            rfe <- fit.rfe(x = data, y = labels, ...)

            # Sort probes
            final <- rfe$features

            # Use "make.names" key to return to original row.names
            final <- merge(data.frame("new" = final), key)$old

            if(length(final) < 2) stop("Uh oh! fsPathClassRFE did not find enough features!")
            array <- new("ExprsBinary",
                         exprs = object@exprs[final, ],
                         annot = object@annot,
                         preFilter = append(object@preFilter, list(final)),
                         reductionModel = append(object@reductionModel, list(NA))
            )

            return(array)
          }
)

setMethod("fsEbayes", "ExprsBinary",
          function(object, probes, ...){ #args to ebayes

            require(limma)

            # Convert 'numeric' probe argument to 'character' probe vector
            if(class(probes) == "numeric"){

              if(probes == 0) probes <- nrow(object@exprs)
              probes <- rownames(object@exprs[1:probes, ])
            }

            # Build data using supplied 'character' probe vector
            if(class(probes) == "character"){

              data <- t(object@exprs[probes, ])
            }

            # Set up and perform eBayes
            design <- as.matrix(ifelse(object@annot$defineCase == "Case", 1, 0))
            colnames(design) <- "CaseVCont"
            fit <- lmFit(object@exprs[probes, ], design)
            ebaye <- ebayes(fit, ...)

            # Sort probes
            vals <- data.frame("probe" = rownames(ebaye$p.value), "p.value" = ebaye$p.value[, 1])
            final <- as.character(vals[order(vals$p.value), "probe"])

            array <- new("ExprsBinary",
                         exprs = object@exprs[final,],
                         annot = object@annot,
                         preFilter = append(object@preFilter, list(final)),
                         reductionModel = append(object@reductionModel, list(NA))
            )

            return(array)
          }
)

# NOTE: mRMR.classic crashes when supplied a very large 'feature_count' argument
setMethod("fsMrmre", "ExprsBinary",
          function(object, probes, ...){ # args to mRMR.classic

            require(mRMRe)

            args <- as.list(substitute(list(...)))[-1]

            if(!"target_indices" %in% names(args)){

              cat("Setting 'target_indices' to 1 (default behavior, override explicitly)...\n")
              args <- append(args, list("target_indices" = 1))
            }

            if(!"feature_count" %in% names(args)){

              cat("Setting 'feature_count' to 64 (default behavior, override explicitly)...\n")
              args <- append(args, list("feature_count" = 64))
            }

            # Convert 'numeric' probe argument to 'character' probe vector
            if(class(probes) == "numeric"){

              if(probes == 0) probes <- nrow(object@exprs)
              probes <- rownames(object@exprs[1:probes, ])
              data <- t(object@exprs[probes, ])
            }

            # Build data using supplied 'character' probe vector
            if(class(probes) == "character"){

              data <- t(object@exprs[probes, ])
            }

            # Set up "make.names" key for improper @exprs row.names
            key <- data.frame("old" = colnames(data), "new" = make.names(colnames(data)))

            # Set up and perform mRMR
            labels <- as.numeric(object@annot$defineCase == "Case")
            mRMRdata <- mRMR.data(data = data.frame(labels, data))
            args <- append(list("data" = mRMRdata), args)
            mRMRout <- do.call(mRMR.classic, args)

            # Sort probes
            final <- as.vector(apply(solutions(mRMRout)[[1]], 2, function(x, y) { return(y[x]) }, y = mRMRe::featureNames(mRMRdata)))

            # Use "make.names" key to return to original row.names
            final <- merge(data.frame("new" = final), key)$old

            array <- new("ExprsBinary",
                         exprs = object@exprs[final,],
                         annot = object@annot,
                         preFilter = append(object@preFilter, list(final)),
                         reductionModel = append(object@reductionModel, list(NA))
            )

            return(array)

          }
)
