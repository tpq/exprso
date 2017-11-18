#' Workhorse for fs Methods
#'
#' Used as a back-end wrapper for creating new fs methods.
#'
#' If the uniqueFx returns a character vector, it is assumed
#'  that the fs method is for feature selection only. If the
#'  uniqueFx returns a list, it is assumed that the fs method
#'  is a reduction model method only.
#'
#' @param object Specifies the \code{ExprsArray} object to undergo feature selection.
#' @param top A numeric scalar or character vector. A numeric scalar indicates
#'  the number of top features that should undergo feature selection. A character vector
#'  indicates specifically which features by name should undergo feature selection.
#'  Set \code{top = 0} to include all features. A numeric vector can also be used
#'  to indicate specific features by location, similar to a character vector.
#' @param uniqueFx A function call unique to the method.
#' @param ... Arguments passed to the detailed function.
#' @return Returns an \code{ExprsArray} object.
#' @export
fs. <- function(object, top, uniqueFx, ...){

  # Convert top input to explicit feature reference
  if(class(top) == "numeric"){
    if(length(top) == 1){
      if(top > nrow(object@exprs)) top <- 0
      if(top == 0) top <- nrow(object@exprs)
      top <- rownames(object@exprs[1:top, ])
    }else{
      top <- rownames(object@exprs[top, ])
    }
  }

  # Build data from top subset for uniqueFx
  data <- t(object@exprs[top, ])
  outcome <- object@annot$defineCase
  final <- do.call("uniqueFx", list(data, outcome, top, ...))

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
fsSample <- function(object, top = 0, ...){ # args to sample

  classCheck(object, "ExprsArray",
             "This function is applied to the results of ?exprso.")

  fs.(object, top,
      uniqueFx = function(data, outcome, top, ...){

        sample(top, ...)
      }, ...)
}

#' Null Feature Selection
#'
#' \code{fsNULL} selects features by interpreting the \code{top} argument.
#'
#' @inheritParams fs.
#' @return Returns an \code{ExprsArray} object.
#' @export
fsNULL <- function(object, top = 0, ...){ # args to NULL

  classCheck(object, "ExprsArray",
             "This function is applied to the results of ?exprso.")

  fs.(object, top,
      uniqueFx = function(data, outcome, top, ...){

        top
      }, ...)
}

#' Select Features by Explicit Reference
#'
#' \code{fsInclude} selects features passed to the \code{include} argument.
#'
#' @inheritParams fs.
#' @param include A character vector. The names of features to rank above all others.
#'  This preserves the feature order otherwise. Argument for \code{fsInclude} only.
#' @return Returns an \code{ExprsArray} object.
#' @export
fsInclude <- function(object, top = 0, include){

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
        c(colnames(data)[index], colnames(data)[!index])
      }, include)
}

#' Select Features by ANOVA
#'
#' \code{fsANOVA} selects features using the \code{aov} function.
#'
#' @inheritParams fs.
#' @return Returns an \code{ExprsArray} object.
#' @export
fsANOVA <- function(object, top = 0, ...){ # args to aov

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
      }, ...)
}

#' Select Features by Statistical Testing
#'
#' \code{fsStats} selects features using the \code{t.test} or \code{ks.test}
#'  function (toggled by the \code{how} argument).
#'
#' @inheritParams fs.
#' @param how A character string. Toggles between the "t.test" and
#'  "ks.test" method. Argument for \code{fsStats} only.
#' @return Returns an \code{ExprsArray} object.
#' @export
fsStats <- function(object, top = 0, how = c("t.test", "ks.test"), ...){ # args to t.test or ks.test

  classCheck(object, "ExprsBinary",
             "This feature selection method only works for binary classification tasks.")

  if(how[1] == "t.test"){

    fs.(object, top,
        uniqueFx = function(data, outcome, top, ...){

          # Prepare data for statistical tests
          cases <- outcome %in% "Case"
          conts <- outcome %in% "Control"
          p <- vector("numeric", length(top))

          for(i in 1:length(top)){
            tryCatch({
              p[i] <- t.test(data[cases, top[i]],
                             data[conts, top[i]], ...)$p.value
            }, error = function(e){
              cat("fsStats failed for feature: ", top[i], ". Setting p(x)=1...\n")
              p[i] <- 1
            })
          }
          top[order(p)]
        }, ...)

  }else if(how[1] == "ks.test"){

    fs.(object, top,
        uniqueFx = function(data, outcome, top, ...){

          # Prepare data for statistical tests
          cases <- outcome %in% "Case"
          conts <- outcome %in% "Control"
          p <- vector("numeric", length(top))

          for(i in 1:length(top)){
            tryCatch({
              p[i] <- ks.test(data[cases, top[i]],
                              data[conts, top[i]], ...)$p.value
            }, error = function(e){
              cat("fsStats failed for feature: ", top[i], ". Setting p(x)=1...\n")
              p[i] <- 1
            })
          }
          top[order(p)]
        }, ...)

  }else{

    stop("Uh oh! Provided 'how' argument not recognized.")
  }
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
fsPrcomp <- function(object, top = 0, ...){ # args to prcomp

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
      }, ...)
}

#' Select Features by Recursive Feature Elimination
#'
#' \code{fsPathClassRFE} selects features using the \code{fit.rfe} function
#'  from the \code{pathClass} package.
#'
#' @inheritParams fs.
#' @return Returns an \code{ExprsArray} object.
#' @export
fsPathClassRFE <- function(object, top = 0, ...){ # args to fit.rfe

  packageCheck("pathClass")
  classCheck(object, "ExprsBinary",
             "This feature selection method only works for binary classification tasks.")

  fs.(object, top,
      uniqueFx = function(data, outcome, top, ...){

        # Set up "make.names" key for improper @exprs row.names
        labels <- factor(outcome, levels = c("Control", "Case"))
        key <- data.frame("old" = colnames(data), "new" = make.names(colnames(data)))

        # NOTE: RFE as assembled by pathClass is via a linear kernel only
        # NOTE: By default, fit.rfe iterates through C = 10^c(-3:3)
        # Run fit.rfe()
        rfe <- pathClass::fit.rfe(x = data, y = labels, ...)

        # Use "make.names" key to return to original row.names
        final <- merge(data.frame("new" = rfe$features), key, sort = FALSE)$old
        if(length(final) < 2) stop("Uh oh! fsPathClassRFE did not find enough features!")
        as.character(final)
      }, ...)
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
fsEbayes <- function(object, top = 0, ...){ # args to ebayes

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
      }, ...)
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
fsEdger <- function(object, top = 0, ...){ # args to exactTest

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
      }, ...)
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
fsMrmre <- function(object, top = 0, ...){ # args to mRMR.classic

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
      }, ...)
}

#' Select Features by Differential Proportionality Analysis
#'
#' \code{fsPropd} selects features using the \code{propd} function
#'  from the \code{propr} package.
#'
#' @inheritParams fs.
#' @return Returns an \code{ExprsArray} object.
#' @export
fsPropd <- function(object, top = 0){

  packageCheck("propr")
  classCheck(object, "ExprsBinary",
             "This feature selection method only works for binary classification tasks.")

  fs.(object, top,
      uniqueFx = function(data, outcome, top){

        # Order pairs by theta
        pd <- propr::propd(data, outcome)
        pd@theta <- pd@theta[order(pd@theta$theta),]

        # Index features by when they first appear
        nrows <- nrow(pd@theta)
        index <- floor(seq(1, nrows+.5, .5))
        odds <- as.logical(1:(nrows*2) %% 2)
        index[odds] <- index[odds] + nrows
        join <- c(pd@theta$Partner, pd@theta$Pair)
        join <- join[index]

        # Rank features by first appearance
        rankedfeats <- unique(join)
        top[rankedfeats]
      })
}
