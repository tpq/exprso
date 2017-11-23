#' Replicate Data Process History
#'
#' \code{modHistory} replicates the \code{fs} history of a reference object.
#'  Used by \code{predict} to prepare validation set for model deployment.
#'
#' @param object An \code{ExprsArray} object. The object that should undergo a
#'  replication of some feature selection and dimension reduction history.
#' @param reference An \code{ExprsArray} or \code{ExprsModel} object. The object
#'  containing the history to use as a template.
#' @return A pre-processed \code{ExprsArray} object.
#' @export
modHistory <- function(object, reference){

  classCheck(object, "ExprsArray",
             "This function is applied to the results of ?exprso.")

  if(!is.null(object@preFilter)){
    if(length(reference@preFilter) <= length(object@preFilter)){
      stop("The object has more history than the provided reference.") }
    if(!identical(object@preFilter, reference@preFilter[1:length(object@preFilter)])){
      stop("The object history does not match the reference history.") }
  }

  # Duplicate starting with first non-overlapping history
  index <- length(object@reductionModel) + 1
  for(i in index:length(reference@reductionModel)){

    if(any(is.na(reference@reductionModel[[i]]))){ # apply fs only
      exprs.i <- object@exprs[reference@preFilter[[i]], , drop = FALSE]
    }else{ # apply fs and dimension reduction
      data <- data.frame(t(object@exprs[reference@preFilter[[i]], , drop = FALSE]))
      if("prcomp" %in% class(reference@reductionModel[[i]])){
        exprs.i <- t(predict(reference@reductionModel[[i]], data))
      }else{ stop("Reduction model class not recognized.") }
    }

    object <- new(class(object), exprs = exprs.i, annot = object@annot,
                  preFilter = append(object@preFilter, list(reference@preFilter[[i]])),
                  reductionModel = append(object@reductionModel,
                                          list(reference@reductionModel[[i]])))
  }

  return(object)
}

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
#' @param object An \code{ExprsArray} object to undergo pre-processing.
#' @param threshold A numeric scalar. The value below which to assign this value.
#' @param maximum A numeric scalar. The value above which to assign this value.
#' @param beta1 A numeric scalar. The \code{max - min} range above which to
#'  include the feature. Inclusive with \code{beta2}.
#' @param beta2 A numeric scalar. The \code{max / min} ratio above which to
#'  include the feature. Inclusive with \code{beta1}.
#' @return A pre-processed \code{ExprsArray} object.
#' @export
modFilter <- function(object, threshold, maximum, beta1, beta2){

  classCheck(object, "ExprsArray",
             "This function is applied to the results of ?exprso.")

  # Set minimum and maximum threshold
  object@exprs[object@exprs < threshold] <- threshold
  object@exprs[object@exprs > maximum] <- maximum

  # Calculate RANGE and MAX:MIN
  index.beta1 <- apply(object@exprs, MARGIN = 1, function(gene) max(gene) - min(gene))
  index.beta2 <- apply(object@exprs, MARGIN = 1, function(gene) max(gene) / min(gene))

  # INCLUDE those with range GREATER THAN beta1 AND beta2
  object@exprs <- object@exprs[index.beta1 > beta1 & index.beta2 > beta2, ]

  return(object)
}

#' Log Transform Data
#'
#' \code{modTransform} log (base 2) transforms feature data.
#'
#' @inheritParams modFilter
#' @return A pre-processed \code{ExprsArray} object.
#' @export
modTransform <- function(object){

  classCheck(object, "ExprsArray",
             "This function is applied to the results of ?exprso.")

  object@exprs <- log(object@exprs, base = 2)
  return(object)
}

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
#' @export
modNormalize <- function(object, MARGIN = c(1, 2)){

  classCheck(object, "ExprsArray",
             "This function is applied to the results of ?exprso.")

  if(2 %in% MARGIN) object@exprs <- apply(object@exprs, MARGIN = 2,
                                          function(samp) (samp - mean(samp)) / sd(samp))
  if(1 %in% MARGIN) object@exprs <- t(apply(object@exprs, MARGIN = 1,
                                            function(gene) (gene - mean(gene)) / sd(gene)))
  return(object)
}

#' Normalize Data
#'
#' \code{modTMM} normalizes feature data.
#'
#' This method normalizes data using the \code{calcNormFactors} function
#'  from the \code{edgeR} package. It returns the original counts
#'  multiplied by the effective library size factors.
#'
#' @inheritParams modFilter
#' @param method A character string. The method used by \code{calcNormFactors}.
#'  Defaults to the "TMM" method.
#' @return A pre-processed \code{ExprsArray} object.
#' @export
modTMM <- function(object, method = "TMM"){

  packageCheck("edgeR")
  classCheck(object, "ExprsArray",
             "This function is applied to the results of ?exprso.")

  f <- edgeR::calcNormFactors(object@exprs, method = method)
  object@exprs <- t(t(object@exprs) * f)
  return(object)
}

#' Compositionally Constrain Data
#'
#' \code{modAcomp} makes it so that all sample vectors have the same total sum.
#'
#' @inheritParams modFilter
#' @return A pre-processed \code{ExprsArray} object.
#' @export
modAcomp <- function(object){

  classCheck(object, "ExprsArray",
             "This function is applied to the results of ?exprso.")

  object@exprs <- apply(object@exprs, 2, function(x) x / sum(x))
  return(object)
}

#' Log-ratio Transform Data
#'
#' \code{modCLR} applies a centered log-ratio transformation to the data.
#'
#' @inheritParams modFilter
#' @return A pre-processed \code{ExprsArray} object.
#' @export
modCLR <- function(object){

  classCheck(object, "ExprsArray",
             "This function is applied to the results of ?exprso.")

  logX <- log(object@exprs)
  object@exprs <- apply(logX, 2, function(x) x / mean(x))
  return(object)
}
