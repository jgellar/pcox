% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/cf.R
\name{cf}
\alias{cf}
\title{Define concurrent time-varying covariate term in a pcox formula}
\usage{
cf(..., lag = 0, linear = TRUE, tv = FALSE, basistype = c("s", "te",
  "t2"), sind = NULL)
}
\arguments{
\item{...}{a list of variables that are the covariates used in the term, as well
as possibly additional arguments that are passed onto the basis constructor
defined by \code{basistype}. Concurrent time-varying covariates should be
included as an \eqn{N x J} matrix, where \eqn{N} is the number of subjects
and \eqn{J} is the number of time points. The right-most columns of most rows
will have \code{NA} values if they are after the event/censoring time.}

\item{lag}{optional time lag for the concurrent effect. Term will be processed
as \eqn{X_i(t-lag)}. Defaults to 0 (no lag).}

\item{linear}{if \code{FALSE}, covariates are included as nonlinear (smooth)
effects, otherwise as a linear effect}

\item{tv}{if \code{TRUE}, makes the effect time-varying}

\item{basistype}{character string that specifies the basis constructor
function (from the \code{mgcv} package) that is used to define a smooth
term. For concurrent terms this is only relevant if the term is
time-varying and nonlinear, because smooths of one variable must use
\code{s}.}

\item{sind}{specifies the time indices for the time-varying covariate. May
  be entered as a vector of length \code{ncol(X)}, or a matrix of the same
  dimensions as \code{X} (for covariates measured on unequal grids).}
}
\value{
The result of a call to \code{p()}, which will be a list with
  the raw data required for the term, and a function of \eqn{x} and \eqn{t}
  that specifies how to set up the term within \code{coxph()}.
}
\description{
Function used to define a concurrent time-varying covariate term in a
pcox formula. Simply calls \code{\link{p}} with default arguments for
concurrent terms.
}
\details{
A concurrent functional term in a \code{pcox} formula is a term
  that involves \eqn{X_i(t)}, or alternatively \eqn{X_i(t-lag)}. This term
  may have a linear or nonlinear effect of \eqn{X_i(t)} and may be time-fixed
  or time-varying (e.g., \eqn{\beta(t)X_i(t)}).
}
\author{
Jonathan Gellar <jgellar1@jhu.edu>
}
\seealso{
\code{\link{p}}
}

