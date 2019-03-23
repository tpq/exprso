#' Workhorse for fs Methods
#'
#' Used as a back-end wrapper for creating new fs methods.
#'
#' If the uniqueFx returns a character vector, it is assumed
#'  that the fs method is for feature selection only. If the
#'  uniqueFx returns a list, it is assumed that the fs method
#'  is a reduction model method only.
#'
#' @param object An \code{ExprsArray} object to undergo feature selection.
#' @param top A numeric scalar or character vector. A numeric scalar indicates
#'  the number of top features that should undergo feature selection. A character vector
#'  indicates specifically which features by name should undergo feature selection.
#'  Set \code{top = 0} to include all features. A numeric vector can also be used
#'  to indicate specific features by location, similar to a character vector.
#' @param uniqueFx A function call unique to the method.
#' @param keep A numeric scalar. Specifies the number of top features that should get
#'  returned by the feature selection method. Use of \code{keep} is generally not
#'  recommended, but can speed up analyses of large data.
#' @param ... Arguments passed to the detailed function.
#' @return Returns an \code{ExprsArray} object.
#' @export
fs. <- function(object, top, uniqueFx, keep, ...){

  # Convert top input to explicit feature reference
  if(class(top) == "numeric"){
    if(length(top) == 1){
      if(top > nrow(object@exprs)) top <- 0
      if(top == 0){
        top <- rownames(object@exprs) # keep for uniqueFx
        data <- t(object@exprs)
      }else{
        top <- rownames(object@exprs)[1:top]
        data <- t(object@exprs[top, , drop = FALSE])
      }
    }else{
      top <- rownames(object@exprs)[top] # when top is numeric vector
      data <- t(object@exprs[top, , drop = FALSE])
    }
  }else{
    # top <- top # feature names already specified
    data <- t(object@exprs[top, , drop = FALSE])
  }

  # Build data from top subset for uniqueFx
  outcome <- object@annot$defineCase
  final <- do.call("uniqueFx", list(data, outcome, top, ...))
  if(keep != 0) final <- final[1:keep]

  # Append uniqueFx results to object
  if(class(final) == "character"){ # fill @preFilter slot
    array <- new(class(object), exprs = object@exprs[final,], annot = object@annot,
                 preFilter = append(object@preFilter, list(final)),
                 reductionModel = append(object@reductionModel, list(NA)))
  }else if(class(final) == "list"){ # fill @reductionModel slot
    array <- new(class(object), exprs = final[[1]], annot = object@annot,
                 preFilter = append(object@preFilter, list(top)),
                 reductionModel = append(object@reductionModel, list(final[[2]])))
  }else{ stop("Uh oh! DEBUG ERROR: 002")}

  return(array)
}

#' Select Features by Random Sampling
#'
#' \code{fsSample} selects features using the \code{sample} function.
#'
#' @inheritParams fs.
#' @return Returns an \code{ExprsArray} object.
#' @export
fsSample <- function(object, top = 0, keep = 0, ...){ # args to sample

  classCheck(object, "ExprsArray",
             "This function is applied to the results of ?exprso.")

  fs.(object, top,
      uniqueFx = function(data, outcome, top, ...){

        sample(top, ...)
      }, keep, ...)
}

#' Null Feature Selection
#'
#' \code{fsNULL} selects features by passing along the \code{top} argument.
#'
#' @inheritParams fs.
#' @return Returns an \code{ExprsArray} object.
#' @export
fsNULL <- function(object, top = 0, keep = 0, ...){ # args to NULL

  classCheck(object, "ExprsArray",
             "This function is applied to the results of ?exprso.")

  fs.(object, top,
      uniqueFx = function(data, outcome, top, ...){

        top
      }, keep, ...)
}

#' Select Features by Explicit Reference
#'
#' \code{fsInclude} selects features passed to the \code{include} argument.
#'  Ranks features by the provided order.
#'
#' @inheritParams fs.
#' @param include A character vector. The names of features to rank above all others.
#'  This preserves the feature order otherwise. Argument for \code{fsInclude} only.
#' @return Returns an \code{ExprsArray} object.
#' @export
fsInclude <- function(object, top = 0, keep = 0, include){

  classCheck(object, "ExprsArray",
             "This function is applied to the results of ?exprso.")

  fs.(object, top,
      uniqueFx = function(data, outcome, top, include){

        if(class(include) != "character"){
          stop("Uh oh! 'include' requires features specified by name.")
        }

        if(!all(include %in% colnames(data))){
          stop("Uh oh! Not all 'include' found in data.")
        }

        index <- colnames(data) %in% include
        c(include, colnames(data)[!index])
      }, keep, include)
}

#' Select Features by ANOVA
#'
#' \code{fsANOVA} selects features using the \code{aov} function.
#'  Note that ANOVA assumes equal variances.
#'
#' @inheritParams fs.
#' @return Returns an \code{ExprsArray} object.
#' @export
fsANOVA <- function(object, top = 0, keep = 0, ...){ # args to aov

  classCheck(object, c("ExprsBinary", "ExprsMulti"),
             "This feature selection method only works for classification tasks.")

  fs.(object, top,
      uniqueFx = function(data, outcome, top, ...){

        # Perform an ANOVA for each feature in data
        df <- data.frame(data, "label" = outcome)
        p <- vector("numeric", length(top))
        for(i in 1:length(top)){

          formula <- stats::as.formula(paste(top[i], "~", "label"))
          fit <- stats::aov(formula, data = df, ...)
          p[i] <- summary(fit)[[1]][1, "Pr(>F)"]
        }

        top[order(p)]
      }, keep, ...)
}

#' Select Features by Statistical Testing
#'
#' \code{fsStats} selects features using a base R statistics
#'  function (toggled by the \code{how} argument).
#'
#' @inheritParams fs.
#' @param how A character string. Toggles between the "t.test", "ks.test",
#'  "wilcox.test", and "var.test" methods.
#' @return Returns an \code{ExprsArray} object.
#' @export
fsStats <- function(object, top = 0, keep = 0,
                    how = c("t.test", "ks.test", "wilcox.test", "var.test"), ...){ # args to base R function

  classCheck(object, "ExprsBinary",
             "This feature selection method only works for binary classification tasks.")

  if(how[1] == "t.test"){

    fs.(object, top,
        uniqueFx = function(data, outcome, top, ...){

          dca <- data[outcome == "Case", ]
          dco <- data[outcome == "Control", ]
          p <- vector("numeric", length(top))
          for(i in 1:ncol(data)){
            tryCatch({
              p[i] <- t.test(dca[, i], dco[, i], ...)$p.value
            }, error = function(e){
              cat("fsStats failed for feature: ", top[i], ". Setting p(x)=1...\n")
              p[i] <- 1
            })
          }
          top[order(p)]
        }, keep, ...)

  }else if(how[1] == "ks.test"){

    fs.(object, top,
        uniqueFx = function(data, outcome, top, ...){

          dca <- data[outcome == "Case", ]
          dco <- data[outcome == "Control", ]
          p <- vector("numeric", length(top))
          for(i in 1:ncol(data)){
            tryCatch({
              p[i] <- ks.test(dca[, i], dco[, i], ...)$p.value
            }, error = function(e){
              cat("fsStats failed for feature: ", top[i], ". Setting p(x)=1...\n")
              p[i] <- 1
            })
          }
          top[order(p)]
        }, keep, ...)

  }else if(how[1] == "wilcox.test"){

    fs.(object, top,
        uniqueFx = function(data, outcome, top, ...){

          dca <- data[outcome == "Case", ]
          dco <- data[outcome == "Control", ]
          p <- vector("numeric", length(top))
          for(i in 1:ncol(data)){
            tryCatch({
              p[i] <- wilcox.test(dca[, i], dco[, i], ...)$p.value
            }, error = function(e){
              cat("fsStats failed for feature: ", top[i], ". Setting p(x)=1...\n")
              p[i] <- 1
            })
          }
          top[order(p)]
        }, keep, ...)

  }else if(how[1] == "var.test"){

    fs.(object, top,
        uniqueFx = function(data, outcome, top, ...){

          dca <- data[outcome == "Case", ]
          dco <- data[outcome == "Control", ]
          p <- vector("numeric", length(top))
          for(i in 1:ncol(data)){
            tryCatch({
              p[i] <- var.test(dca[, i], dco[, i], ...)$p.value
            }, error = function(e){
              cat("fsStats failed for feature: ", top[i], ". Setting p(x)=1...\n")
              p[i] <- 1
            })
          }
          top[order(p)]
        }, keep, ...)

  }else{

    stop("Uh oh! Provided 'how' argument not recognized.")
  }
}

#' Select Features by Correlation
#'
#' \code{fsCor} selects features using the \code{cor} function.
#'  Ranks features by absolute value of correlation.
#'
#' @inheritParams fs.
#' @return Returns an \code{ExprsArray} object.
#' @export
fsCor <- function(object, top = 0, keep = 0, ...){ # args to cor

  classCheck(object, "RegrsArray",
             "This feature selection method only works for continuous outcome tasks.")

  fs.(object, top,
      uniqueFx = function(data, outcome, top, ...){

        r <- vector("numeric", length(top))
        for(i in 1:ncol(data)){
          r[i] <- cor(data[,i], outcome)
        }
        top[rev(order(r))]
      }, keep, ...)
}

#' Reduce Dimensions by PCA
#'
#' \code{fsPrcomp} reduces dimensions using the \code{prcomp} function.
#'  The reduction model is saved and deployed automatically on any new
#'  data during model validation.
#'
#' @inheritParams fs.
#' @return Returns an \code{ExprsArray} object.
#' @export
fsPrcomp <- function(object, top = 0, keep = 0, ...){ # args to prcomp

  classCheck(object, "ExprsArray",
             "This function is applied to the results of ?exprso.")

  fs.(object, top,
      uniqueFx = function(data, outcome, top, ...){

        # ATTENTION: We want dependent variables as columns
        # NOTE: as.data.frame will not rename columns
        reductionModel <- prcomp(data.frame(data), ...)

        # ATTENTION: The value of predict(reductionModel, data) equals $x
        # This information will automatically distill the data
        #  when calling modHistory
        list(t(reductionModel$x),
             reductionModel)
      }, keep, ...)
}

#' Reduce Dimensions by PCA
#'
#' \code{fsPCA} reduces dimensions using the \code{prcomp} function.
#'  The reduction model is saved and deployed automatically on any new
#'  data during model validation.
#'
#' @inheritParams fs.
#' @return Returns an \code{ExprsArray} object.
#' @export
fsPCA <- function(object, top = 0, keep = 0, ...){ # args to prcomp

  fsPrcomp(object, top = top, keep = keep, ...)
}

#' Reduce Dimensions by RDA
#'
#' \code{fsRDA} reduces dimensions using the \code{rda} function
#'  from the \code{vegan} package. The reduction model is saved and
#'  deployed automatically on any new data during model validation.
#'
#' @inheritParams fs.
#' @return Returns an \code{ExprsArray} object.
#' @export
fsRDA <- function(object, top = 0, keep = 0, ...){ # args to rda

  packageCheck("vegan")
  classCheck(object, "ExprsArray",
             "This function is applied to the results of ?exprso.")

  fs.(object, top,
      uniqueFx = function(data, outcome, top, ...){

        # ATTENTION: We want dependent variables as columns
        # NOTE: as.data.frame will not rename columns
        reductionModel <- vegan::rda(data.frame(data), ...)

        # ATTENTION: predict(reductionModel, data, type = "wa") equals $CA$u
        list(t(reductionModel$CA$u),
             reductionModel)
      }, keep, ...)
}

#' Select Features by Moderated t-test
#'
#' \code{fsEbayes} selects features using the \code{lmFit} and
#'  \code{eBayes} functions from the \code{limma} package. Features
#'  ranked by the \code{topTableF} function.
#'
#' @inheritParams fs.
#' @return Returns an \code{ExprsArray} object.
#' @export
fsEbayes <- function(object, top = 0, keep = 0, ...){ # args to ebayes

  packageCheck("limma")
  classCheck(object, c("ExprsBinary", "ExprsMulti"),
             "This feature selection method only works for classification tasks.")

  fs.(object, top,
      uniqueFx = function(data, outcome, top, ...){

        design <- stats::model.matrix(~0 + outcome)
        fit <- limma::lmFit(t(data), design)
        ebaye <- limma::eBayes(fit)
        tt <- limma::topTableF(ebaye, number = ncol(data))
        rownames(tt)
      }, keep, ...)
}

#' Selects Features by Exact Test
#'
#' \code{fsEdger} selects features using the \code{exactTest} function
#'  from the \code{edgeR} package. This function does not normalize the data,
#'  but does estimate dispersion using the \code{estimateCommonDisp}
#'  and \code{estimateTagwiseDisp} functions.
#'
#' The user can normalize the data before feature selection using the
#'  \code{modTMM} function. Note that applying \code{edgeR} to already normalized
#'  counts differs slightly from applying \code{edgeR} with normalization.
#'
#' @inheritParams fs.
#' @return Returns an \code{ExprsArray} object.
#' @export
fsEdger <- function(object, top = 0, keep = 0, ...){ # args to exactTest

  packageCheck("edgeR")
  classCheck(object, "ExprsBinary",
             "This feature selection method only works for binary classification tasks.")

  fs.(object, top,
      uniqueFx = function(data, outcome, top, ...){

        y <- edgeR::DGEList(counts = t(data), group = outcome)
        y <- edgeR::estimateCommonDisp(y)
        y <- edgeR::estimateTagwiseDisp(y)
        et <- edgeR::exactTest(y)
        tt <- as.data.frame(edgeR::topTags(et, n = nrow(et)))
        rownames(tt)
      }, keep, ...)
}

#' Select Features by mRMR
#'
#' \code{fsMrmre} selects features using the \code{mRMR.classic} function
#'  from the \code{mRMRe} package.
#'
#' Note that \code{fsMrmre} crashes when supplied a very large
#'  \code{feature_count} owing to its \code{mRMRe} implementation.
#'
#' @inheritParams fs.
#' @return Returns an \code{ExprsArray} object.
#' @export
fsMrmre <- function(object, top = 0, keep = 0, ...){ # args to mRMR.classic

  packageCheck("mRMRe")
  classCheck(object, "ExprsBinary",
             "This feature selection method only works for binary classification tasks.")

  fs.(object, top,
      uniqueFx = function(data, outcome, top, ...){

        args <- getArgs(...)
        args <- defaultArg("target_indices", 1, args)
        args <- defaultArg("feature_count", 64, args)

        # Set up "make.names" key for improper @exprs row.names
        key <- data.frame("old" = colnames(data), "new" = make.names(colnames(data)))

        labels <- as.numeric(outcome == "Case")
        mRMRdata <- mRMRe::mRMR.data(data = data.frame(labels, data))
        args <- append(list("data" = mRMRdata), args)
        mRMRout <- do.call(mRMRe::mRMR.classic, args)

        # Sort features
        final <- as.vector(
          apply(mRMRe::solutions(mRMRout)[[1]], 2, function(x, y){ return(y[x])},
                y = mRMRe::featureNames(mRMRdata))
        )

        # Use "make.names" key to return to original row.names
        final <- merge(data.frame("new" = final), key, sort = FALSE)$old
        as.character(final)
      }, keep, ...)
}

#' Select Features by Rank Product Analysis
#'
#' \code{fsRankProd} selects features using the \code{RankProducts} function
#'  from the \code{RankProd} package.
#'
#' @inheritParams fs.
#' @return Returns an \code{ExprsArray} object.
#' @export
fsRankProd <- function(object, top = 0, keep = 0, ...){ # args to RankProducts

  packageCheck("RankProd")
  classCheck(object, "ExprsBinary",
             "This feature selection method only works for binary classification tasks.")

  fs.(object, top,
      uniqueFx = function(data, outcome, top, ...){

        args <- getArgs(...)
        args <- defaultArg("logged", FALSE, args)

        cl <- ifelse(outcome == "Control", 0, 1)
        args <- append(list("data" = t(data), "cl" = cl), args)
        sink(tempfile())
        final <- do.call(RankProd::RankProducts, args)
        sink()

        # Get two-tailed p-value
        p <- apply(final$pval, 1, min)
        top[order(p)]
      }, keep, ...)
}

#' Convert Features into Balances
#'
#' \code{fsBalance} converts features into balances.
#'
#' @inheritParams fs.
#' @param sbp.how A character string. The method used to build
#'  the serial binary partition matrix of balances. Any
#'  \code{balance::sbp.from*} function will work.
#' @param ternary A boolean. Toggles whether to return balances
#'  representing three components. Argument passed to
#'  \code{balance::sbp.subset}. Set \code{ternary = FALSE} and
#'  \code{ratios = FALSE} to skip subset.
#' @param ratios A boolean. Toggles whether to return balances
#'  representing two components. Argument passed to
#'  \code{balance::sbp.subset}. Set \code{ternary = FALSE} and
#'  \code{ratios = FALSE} to skip subset.
#' @return Returns an \code{ExprsArray} object.
#' @export
fsBalance <- function(object, top = 0, keep = 0, sbp.how = "sbp.fromPBA",
                      ternary = FALSE, ratios = FALSE, ...){ # args to sbp.how

  packageCheck("balance")
  classCheck(object, "ExprsArray",
             "This function is applied to the results of ?exprso.")

  fs.(object, top,
      uniqueFx = function(data, outcome, top, sbp.how, ternary, ratios, ...){

        args <- getArgs(...)
        args <- append(list("x" = data), args)

        if(sbp.how == "sbp.fromPropd"){
          message("Setting group (forced behavior, cannot override)...\n")
          args <- suppressMessages(forceArg("group", outcome, args))
        }

        sbp <- do.call(get(sbp.how, asNamespace("balance")), args)
        sbp <- balance::sbp.subset(sbp, ternary, ratios)
        balances <- t(balance::balance.fromSBP(data, sbp))
        colnames(balances) <- rownames(data)
        class(sbp) <- "SBP"

        list(balances,
             sbp)
      }, keep, sbp.how, ternary, ratios, ...)
}
