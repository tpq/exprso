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
#' @param uniqueFx A function call unique to the method. See Details.
#' @param ... Arguments passed to the detailed function.
#' @return Returns an \code{ExprsArray} object.
#' @export
fs. <- function(object, top, uniqueFx, ...){

  # If a reduction model has never been used, there is no reason to store old history
  # (helps to reduce RAM overhead and to improve run-time)
  if(!is.null(object@reductionModel)){
    if(all(unlist(lapply(object@reductionModel, is.na)))){
      object@preFilter <- lapply(object@preFilter, function(x) 0)
    }
  }

  # Convert top input to explicit feature reference
  if(class(top) == "numeric"){
    if(length(top) == 1){
      if(top > nrow(object@exprs)) top <- 0
      if(top == 0){
        topChar <- rownames(object@exprs) # keep for uniqueFx
      }else{
        topChar <- rownames(object@exprs)[1:top] # when top is numeric scalar
      }
    }else{
      topChar <- rownames(object@exprs)[top] # when top is numeric vector
    }
  }else{
    topChar <- top # else top is row names
  }

  # Run uniqueFx on top data
  if(!identical(top, 0)){
    x <- t(object@exprs[topChar, , drop = FALSE])
  }else{
    x <- t(object@exprs) # runs faster
  }
  y <- object@annot$defineCase
  final <- do.call("uniqueFx", list(data = x, outcome = y, top = topChar, ...))

  # Append uniqueFx results to object
  if(class(final) == "character"){ # fill @preFilter slot

    output <- object
    output@exprs <- object@exprs[final,]
    output@annot <- object@annot
    output@preFilter <- append(object@preFilter, list(final))
    output@reductionModel <- append(object@reductionModel, list(NA))

  }else if(class(final) == "list"){ # fill @reductionModel slot

    output <- object
    output@exprs <- final[[1]]
    output@annot <- object@annot
    output@preFilter <- append(object@preFilter, list(topChar))
    output@reductionModel <- append(object@reductionModel, list(final[[2]]))

  }else{ stop("Uh oh! DEBUG ERROR: FS1")}

  return(output)
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
#' \code{fsNULL} does not select features. However, it will handle the
#'  \code{top} and \code{keep} arguments if provided.
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
#' \code{fsInclude} moves features named by the \code{include} argument
#'  to the top of the ranked features list. Otherwise, the relative order
#'  of the features does not change.
#'
#' @inheritParams fs.
#' @param include A character vector. The names of features to rank above all others.
#' @return Returns an \code{ExprsArray} object.
#' @export
fsInclude <- function(object, top = 0, include){

  classCheck(object, "ExprsArray",
             "This function is applied to the results of ?exprso.")

  if(!identical(top, 0)){
    stop("fsInclude must have top = 0!")
  }

  fs.(object, top = 0,
      uniqueFx = function(data, outcome, top, include){

        if(class(include) != "character"){
          stop("Uh oh! 'include' requires features specified by name.")
        }

        if(!all(include %in% top)){
          stop("Uh oh! Not all 'include' found in data.")
        }

        index <- top %in% include
        c(include, top[!index])
      }, include)
}

#' Select Features by ANOVA
#'
#' \code{fsANOVA} selects features using the \code{aov} function.
#'  Note that the ANOVA assumes equal variances, so will differ from
#'  the \code{fsStats} t-test in the two-group setting.
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
#' \code{fsStats} selects features using a base R statistics
#'  function (toggled by the \code{how} argument).
#'
#' @inheritParams fs.
#' @param how A character string. Toggles between the "t.test", "ks.test",
#'  "wilcox.test", and "var.test" methods.
#' @return Returns an \code{ExprsArray} object.
#' @export
fsStats <- function(object, top = 0,
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
        }, ...)

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
        }, ...)

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
        }, ...)

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
        }, ...)

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
fsCor <- function(object, top = 0, ...){ # args to cor

  classCheck(object, "RegrsArray",
             "This feature selection method only works for continuous outcome tasks.")

  fs.(object, top,
      uniqueFx = function(data, outcome, top, ...){

        r <- vector("numeric", length(top))
        for(i in 1:ncol(data)){
          r[i] <- cor(data[,i], outcome)
        }
        r <- abs(r)
        top[order(r, decreasing = TRUE)]
      }, ...)
}

#' Reduce Dimensions by PCA
#'
#' \code{fsPrcomp} runs a PCA using the \code{prcomp} function.
#'  The PCA model is saved and deployed automatically by the
#'  \code{predict} method during test set validation.
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

#' Reduce Dimensions by PCA
#'
#' \code{fsPrcomp} runs a PCA using the \code{prcomp} function.
#'  The PCA model is saved and deployed automatically by the
#'  \code{predict} method during test set validation.
#'
#' @inheritParams fs.
#' @return Returns an \code{ExprsArray} object.
#' @export
fsPCA <- function(object, top = 0, ...){ # args to prcomp

  fsPrcomp(object, top = top, ...)
}

#' Reduce Dimensions by RDA
#'
#' \code{fsRDA} runs an RDA using the \code{rda} function
#'  from the \code{vegan} package, partialling out \code{colBy}.
#'  The RDA model is saved and deployed automatically by the
#'  \code{predict} method during test set validation.
#'
#' When \code{colBy} is provided, it serves as the constraining matrix.
#'  However, \code{fsRDA} always returns the unconstrained scores.
#'  As such, \code{fsRDA} effectively partials out the contribution
#'  of \code{colBy} to the training set, then uses this rule
#'  to partial out the contribution to the test set too.
#'
#' @inheritParams fs.
#' @param colBy A character vector. Lists the columns in \code{@@annot}
#'  to use as the constraining matrix. Passed to \code{vegan::rda}.
#'  Optional argument. Skip with \code{colBy = NULL}.
#' @return Returns an \code{ExprsArray} object.
#' @export
fsRDA <- function(object, top = 0, colBy = NULL){ # args to rda

  packageCheck("vegan")
  classCheck(object, "ExprsArray",
             "This function is applied to the results of ?exprso.")

  if(is.null(colBy)){
    Y <- NULL
  }else{
    Y <- object@annot[, colBy, drop = FALSE]
  }

  fs.(object, top,
      uniqueFx = function(data, outcome, top, Y){

        # NOTE: If colBy provided, run constrained RDA
        if(is.null(colBy)){
          reductionModel <- vegan::rda(data.frame(data))
        }else{
          reductionModel <- vegan::rda(data.frame(data), Y = Y)
        }

        # ATTENTION: predict(reductionModel, data, type = "wa") equals $CA$u
        list(t(reductionModel$CA$u),
             reductionModel)
      }, Y)
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

#' Select Features by Rank Product Analysis
#'
#' \code{fsRankProd} selects features using the \code{RankProducts} function
#'  from the \code{RankProd} package.
#'
#' @inheritParams fs.
#' @return Returns an \code{ExprsArray} object.
#' @export
fsRankProd <- function(object, top = 0, ...){ # args to RankProducts

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
      }, ...)
}

#' Convert Features into Balances
#'
#' \code{fsBalance} converts features into balances.
#'  The balance rule is saved and deployed automatically by the
#'  \code{predict} method during test set validation.
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
fsBalance <- function(object, top = 0, sbp.how = "sbp.fromPBA",
                      ternary = FALSE, ratios = FALSE, ...){ # args to sbp.how

  packageCheck("balance")
  classCheck(object, "ExprsArray",
             "This function is applied to the results of ?exprso.")

  fs.(object, top,
      uniqueFx = function(data, outcome, top, sbp.how, ternary, ratios, ...){

        if(any(data <= 0)){
          stop("This feature selection method does not work for non-positive values.")
        }

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
      }, sbp.how, ternary, ratios, ...)
}

#' Reduce Dimensions by Amalgamation
#'
#' \code{fsAmalgam} finds a set of explanatory "amalgams" using
#'  the \code{amalgam::amalgam} function. This function expects
#'  a compositional data set that can be reduced by amalgamation.
#'  The resultant "amalgams" are clr- or slr-transformed.
#'  The amalgamation rule is saved and deployed automatically by the
#'  \code{predict} method during test set validation.
#'
#' @inheritParams fs.
#' @return Returns an \code{ExprsArray} object.
#' @export
fsAmalgam <- function(object, top = 0, ...){ # args to amalgam

  packageCheck("amalgam")
  classCheck(object, "ExprsArray",
             "This function is applied to the results of ?exprso.")

  fs.(object, top,
      uniqueFx = function(data, outcome, top, ...){

        if(any(data < 0)){
          stop("This feature selection method does not work for negative values.")
        }

        args <- getArgs(...)
        args <- append(list("x" = data), args)

        # Run the amalgam genetic algorithm
        sink(tempfile())
        model <- do.call(amalgam::amalgam, args)
        sink()

        if(is.null(model$SLR)){

          # If not using SLRs, then CLR the data
          reducedData <- model$amalgams # = original %*% weights
          reducedData <- t(apply(reducedData, 1, function(x) log(x) - mean(log(x)))) # clr-transform
          class(model) <- "amalg-clr"

        }else{

          # Else use the SLRS directly
          reducedData <- model$SLR
          class(model) <- "amalg-slr"
        }

        # Make sure sample names carry through...
        reducedData <- t(reducedData)
        colnames(reducedData) <- rownames(data)
        rownames(reducedData) <- paste0("z", 1:nrow(reducedData))

        list(reducedData,
             model)
      }, ...)
}

#' Use Annotations as Features
#'
#' \code{fsAnnot} moves annotations named by the \code{colBy} argument
#'  to the top of the ranked features list. Otherwise, the relative order
#'  of the features does not change.
#'
#' @inheritParams fs.
#' @param colBy A character vector. The names of annotations to rank above all others.
#' @return Returns an \code{ExprsArray} object.
#' @export
fsAnnot <- function(object, top = 0, colBy){

  classCheck(object, "ExprsArray",
             "This function is applied to the results of ?exprso.")

  if(!identical(top, 0)){
    stop("fsAnnot must have top = 0!")
  }

  if("defineCase" %in% colBy) stop("You can't use the outcome as a feature!")

  Y <- object@annot[, colBy, drop = FALSE]
  for(col in 1:ncol(Y)){ Y[,col] <- as.numeric(Y[,col]) }

  fs.(object, top = 0,
      uniqueFx = function(data, outcome, top, Y){

        model.output <- cbind(Y, data)
        model.rule <- colBy
        class(model.rule) <- "fsAnnot"

        list(t(model.output),
             model.rule)
      }, Y)
}

#' Reduce Dimensions by Log-Ratio Selection
#'
#' \code{fsPRA} finds the most explanatory pairwise log-ratios
#'  using the variable selection method proposed by Michael Greenacre
#'  in "Variable Selection in Compositional Data Analysis Using
#'  Pairwise Logratios", modified to run faster.
#'
#' @inheritParams fs.
#' @return Returns an \code{ExprsArray} object.
#' @export
fsPRA <- function(object, top = 0, ...){ # args to pra

  packageCheck("propr")
  classCheck(object, "ExprsArray",
             "This function is applied to the results of ?exprso.")

  fs.(object, top,
      uniqueFx = function(data, outcome, top, ...){

        if(any(data <= 0)){
          stop("This feature selection method does not work for non-positive values.")
        }

        args <- as.list(substitute(list(...)))[-1]
        args <- append(list("counts" = data), args)

        # Run pra on data
        sink(tempfile())
        res <- do.call(propr::pra, args)
        sink()

        reducedData <- res$Y # note: log(pair / partner)
        colnames(reducedData) <- paste0(res$best$Pair, ".vs.", res$best$Partner)
        model <- as.matrix(res$best) # as.matrix needed to set class
        class(model) <- "pra"

        list(t(reducedData),
             model)

      }, ...)
}
