% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/6.1-predict.R, R/8.2-ens.R
\docType{methods}
\name{exprso-predict}
\alias{exprso-predict}
\alias{predict,ExprsMachine-method}
\alias{predict,ExprsModule-method}
\alias{predict,RegrsModel-method}
\alias{predict,ExprsEnsemble-method}
\title{Deploy Model}
\usage{
\S4method{predict}{ExprsMachine}(object, array, verbose = TRUE)

\S4method{predict}{ExprsModule}(object, array, verbose = TRUE)

\S4method{predict}{RegrsModel}(object, array, verbose = TRUE)

\S4method{predict}{ExprsEnsemble}(object, array, how = "majority",
  verbose = TRUE)
}
\arguments{
\item{object}{An \code{exprso} model.}

\item{array}{An \code{exprso} object. The test data.}

\item{verbose}{A boolean. Argument passed to \code{calcStats}.}

\item{how}{A string. Describes how the ensemble decides. By default, it
uses "majority" voting. However, the user can select "probability" voting
for binary classifier ensembles.}
}
\value{
Returns an \code{exprso} prediction object.
}
\description{
Deploy a model to predict outcomes from the data.
}
\details{
Models can only get deployed on an object of the type used to build
 the model. This method now supports binary classification,
 multi-class classification, and regression.

For regression ensembles, the average outcome is reported. For multi-class
 classifier ensembles, the majority vote is reported. For binary classifier
 ensembles, the majority vote or probability-weighted vote is reported.
 For probability-weighted voting considers the threshold, the average
 "Case" probability is reported. All ties broken randomly.
}
