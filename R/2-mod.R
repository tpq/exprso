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
  unsortedFeatures <- rownames(object@exprs)
  index <- length(object@reductionModel) + 1
  for(i in index:length(reference@reductionModel)){ # for every fs in reference...

    indexedFeatures <- reference@preFilter[[i]]
    indexedModel <- reference@reductionModel[[i]]
    if(any(is.na(indexedModel))){ # apply fs (via top) only

      if(identical(indexedFeatures, unsortedFeatures)){
        exprs.i <- object@exprs # avoid unnecessary subsets
      }else{
        exprs.i <- object@exprs[indexedFeatures, , drop = FALSE]
      }

    }else{ # apply fs (via top) and dimension reduction

      data <- data.frame(t(object@exprs[indexedFeatures, , drop = FALSE]))
      if("prcomp" %in% class(indexedModel)){
        exprs.i <- t(predict(indexedModel, data))
      }else if("SBP" %in% class(indexedModel)){
        packageCheck("balance")
        exprs.i <- t(balance::balance.fromSBP(data, indexedModel))
        colnames(exprs.i) <- rownames(data)
      }else{
        stop("Reduction model not recognized.")
      }
    }

    object <- new(class(object), exprs = exprs.i, annot = object@annot,
                  preFilter = append(object@preFilter, list(indexedFeatures)),
                  reductionModel = append(object@reductionModel, list(indexedModel)))
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
#' \code{modTransform} log transforms feature data.
#'
#' @inheritParams modFilter
#' @param base A numeric scalar. The base of the logarithm.
#' @return A pre-processed \code{ExprsArray} object.
#' @export
modTransform <- function(object, base = exp(1)){

  classCheck(object, "ExprsArray",
             "This function is applied to the results of ?exprso.")

  object@exprs <- log(object@exprs, base = base)
  return(object)
}

#' Sample Features from Data
#'
#' \code{modSample} samples features from a data set randomly without
#'  replacement. When \code{size = 0}, this is equivalent to
#'  \code{fsSample, top = 0}, but much quicker.
#'
#' @inheritParams modFilter
#' @param size A numeric scalar. The number of randomly sampled features
#'  to include in the pre-processed \code{ExprsArray} object.
#' @return A pre-processed \code{ExprsArray} object.
#' @export
modSample <- function(object, size = 0){

  classCheck(object, "ExprsArray",
             "This function is applied to the results of ?exprso.")

  if(size == 0) size <- nrow(object@exprs)
  keep <- sample(1:nrow(object@exprs), size = size)
  object@exprs <- object@exprs[keep,,drop = FALSE]
  return(object)
}

#' Permute Features in Data
#'
#' \code{modPermute} randomly samples each feature in the data
#'  without replacement. This method helps establish a null
#'  model for the purpose of testing the significance of
#'  observed prediction error estimates.
#'
#' @inheritParams modFilter
#' @return A pre-processed \code{ExprsArray} object.
#' @export
modPermute <- function(object){

  classCheck(object, "ExprsArray",
             "This function is applied to the results of ?exprso.")

  oldnames <- colnames(object@exprs)
  object@exprs <- t(apply(object@exprs, 1, sample))
  colnames(object@exprs) <- oldnames
  return(object)
}

#' Select Features from Data
#'
#' \code{modSelect} selects specific features from a data set. Unlike
#'  \code{fsInclude}, this function does not update \code{@@preFilter}
#'  and returns only those features stated by \code{include}.
#'
#' @inheritParams modFilter
#' @param include A character vector. The names of features to include.
#' @return A pre-processed \code{ExprsArray} object.
#' @export
modInclude <- function(object, include = rownames(object@exprs)){

  classCheck(object, "ExprsArray",
             "This function is applied to the results of ?exprso.")

  object@exprs <- object@exprs[include,,drop = FALSE]
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

#' Recast Data as Feature (Log-)Ratios
#'
#' \code{modRatios} recasts a data set with N feature columns as a new
#'  data set with N * (N - 1) / 2 feature (log-)ratio columns.
#'
#' @inheritParams modHistory
#' @param alpha A numeric scalar. This argument guides
#'  a Box-Cox transformation to approximate log-ratios in the
#'  presence of zeros. Skip with \code{NA}.
#' @return A pre-processed \code{ExprsArray} object.
#' @export
modRatios <- function(object, alpha = NA){

  packageCheck("propr")
  classCheck(object, "ExprsArray",
             "This function is applied to the results of ?exprso.")

  object@exprs <- t(propr::ratios(t(object@exprs), alpha = NA))
  return(object)
}

#' Scale Data by Factor Range
#'
#' \code{modScale} scales a data set by making all sample vectors
#'  have the same total sum, then multiplying each sample vector by
#'  a scale factor.
#'
#' If \code{uniform = TRUE}, scale factors are randomly sampled from
#'  the uniform distribution \code{(0, alpha) + 1}. Otherwise, scale
#'  factors are randomly sampled from the normal distribution with
#'  a mean of 0 and standard deviation of \code{alpha}. When using
#'  the normal distribution, these scale factors are transformed by
#'  taking the absolute value then adding one. For this reason,
#'  data are always unscaled when \code{alpha = 0}.
#'
#' @inheritParams modHistory
#' @param alpha An integer. The maximum range of scale factors used
#'  for scaling if \code{uniform = TRUE}. The standard deviation
#'  of the scale factors if \code{uniform = FALSE}. See Details.
#' @param uniform A boolean. Toggles whether to draw scale factors
#'  from a uniform distribution or a normal distribution.
#' @return A pre-processed \code{ExprsArray} object.
#' @export
modScale <- function(object, alpha = 0, uniform = TRUE){

  classCheck(object, "ExprsArray",
             "This function is applied to the results of ?exprso.")

  if(uniform){

    # Draw scale factors from Uniform distribution
    lambda <- stats::runif(ncol(object@exprs), min = 0, max = alpha) + 1

  }else{

    # Draw scale factors from Normal distribution
    lambda <- abs(stats::rnorm(ncol(object@exprs), mean = 0, sd = alpha)) + 1
  }

  # Apply scale to weigh samples
  object <- modAcomp(object)
  newdata <- apply(object@exprs, 1, function(x) x * lambda)
  object@exprs <- t(newdata)
  return(object)
}

#' Skew Data by Factor Range
#'
#' \code{modSkew} skews a data set by making all sample vectors
#'  have the same total sum, introducing a new feature, and then
#'  making all sample vectors again have the same total sum.
#'
#' If \code{uniform = TRUE}, skew factors are randomly sampled from
#'  the uniform distribution \code{(0, alpha) + 1}. Otherwise, skew
#'  factors are randomly sampled from the normal distribution with
#'  a mean of 0 and standard deviation of \code{alpha}. When using
#'  the normal distribution, these skew factors are transformed by
#'  taking the absolute value then adding one. For this reason,
#'  data are always unskewed when \code{alpha = 0}.
#'
#' @inheritParams modHistory
#' @param alpha An integer. The maximum range of skew factors used
#'  for skewing if \code{uniform = TRUE}. The standard deviation
#'  of the skew factors if \code{uniform = FALSE}. See Details.
#' @param uniform A boolean. Toggles whether to draw skew factors
#'  from a uniform distribution or a normal distribution.
#' @return A pre-processed \code{ExprsArray} object.
#' @export
modSkew <- function(object, alpha = 0, uniform = TRUE){

  classCheck(object, "ExprsArray",
             "This function is applied to the results of ?exprso.")

  if(uniform){

    # Draw skew factors from Uniform distribution
    lambda <- stats::runif(ncol(object@exprs), min = 0, max = alpha) + 1

  }else{

    # Draw skew factors from Normal distribution
    lambda <- abs(stats::rnorm(ncol(object@exprs), mean = 0, sd = alpha)) + 1
  }

  # Apply scale to weigh samples
  object <- modAcomp(object)

  # Join gamma and close data
  object@exprs <- rbind(object@exprs, lambda)
  object <- modAcomp(object)
  return(object)
}
