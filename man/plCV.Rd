% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/7.1-plCV.R
\name{plCV}
\alias{plCV}
\title{Perform Simple Cross-Validation}
\usage{
plCV(array, top, how, fold, aucSkip, plCV.acc, ...)
}
\arguments{
\item{array}{Specifies the \code{ExprsArray} object to undergo cross-validation.}

\item{top}{A numeric scalar or character vector. A numeric scalar indicates
the number of top features that should undergo feature selection. A character vector
indicates specifically which features by name should undergo feature selection.
Set \code{top = 0} to include all features. Note that providing a numeric vector
for the \code{top} argument will have \code{plGrid} search across multiple
top features. However, by providing a list of numeric vectors as the \code{top}
argument, the user can force the default handling of numeric vectors.}

\item{how}{A character string. The \code{\link{build}} method to iterate.}

\item{fold}{A numeric scalar. The number of folds for cross-validation.
Set \code{fold = 0} to perform leave-one-out cross-validation. Argument passed
to \code{\link{plCV}}. Set \code{fold = NULL} to skip cross-validation altogether.}

\item{aucSkip}{A logical scalar. Argument passed to \code{\link{calcStats}}.}

\item{plCV.acc}{A string. The performance metric to use. For example,
choose from "acc", "sens", "spec", "prec", "f1", "auc", or any of the
regression specific measures. Argument passed to \code{\link{plCV}}.}

\item{...}{Arguments passed to the \code{how} method. Unlike the \code{build} method,
\code{plGrid} allows multiple parameters for each argument, supplied as a vector.
See Details.}
}
\value{
The average inner-fold cross-validation accuracy.
}
\description{
Calculates v-fold or leave-one-out cross-validation without selecting a new
 set of features with each fold. See Details.
}
\details{
\code{plCV} performs v-fold or leave-one-out cross-validation. The argument
 \code{fold} specifies the number of v-folds to use during cross-validation.
 Set \code{fold = 0} to perform leave-one-out cross-validation.

This type of cross-validation is most appropriate if the data
 has not undergone any prior feature selection. However, it is also useful
 as an unbiased guide to parameter selection within another
 \code{\link{pl}} workflow.

Users should never need to call this function directly. Instead, they
 should use \code{\link{plMonteCarlo}} or \code{\link{plNested}}.
 There, \code{plCV} handles inner-fold cross-validation.
}
