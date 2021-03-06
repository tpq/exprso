% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/5.1-fs.R
\name{fsANOVA}
\alias{fsANOVA}
\title{Select Features by ANOVA}
\usage{
fsANOVA(object, top = 0, ...)
}
\arguments{
\item{object}{An \code{ExprsArray} object to undergo feature selection.}

\item{top}{A numeric scalar or character vector. A numeric scalar indicates
the number of top features that should undergo feature selection. A character vector
indicates specifically which features by name should undergo feature selection.
Set \code{top = 0} to include all features. A numeric vector can also be used
to indicate specific features by location, similar to a character vector.}

\item{...}{Arguments passed to the detailed function.}
}
\value{
Returns an \code{ExprsArray} object.
}
\description{
\code{fsANOVA} selects features using the \code{aov} function.
 Note that the ANOVA assumes equal variances, so will differ from
 the \code{fsStats} t-test in the two-group setting.
}
