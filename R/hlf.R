#' Constuct a Cox regression term as a historical linear function
#' 
#' This function defines a Cox regression term as a historical linear function
#' for inclusion in a \code{pcox} formula. This term treats a densely measured
#' time-varying covariate (i.e., biomarker) as a functional predictor,
#' and models its effect on the hazard at time t with as a smooth
#' funcitonal term. The default is the time-static term
#' \eqn{1\t \int_0^t X_i(s)\beta(s)ds}. A common modification would be
#' to make the term time-varying by specifying it as a\code{\link{tv}}
#' term, resulting in the formula \eqn{1\t \int_0^t X_i(s)\beta(s,t)ds}.
#' The \eqn{1/t} part of the formula may be removed and the integration
#' limits can be changed.
#' 
#' @param X matrix containing the time-varying covariate. Should be
#'    \eqn{N x J}, where \eqn{N} is the number of subjects and \eqn{J}
#'    is the maximum number of time points per subject. Most rows will
#'    have \code{NA} values in the right-most columns, corresponding
#'    to unobserved time points.
#' @param tind matrix containing the time indices of evaluations of
#'    \eqn{X_i(t)}. If a matrix, it must be the same dimensionality
#'    as \code{X}; if a vector, must be of length \code{ncol(X)} and
#'    the same grid is assumed for all subjects; if omitted, the same
#'    equally-spaced grid is assumed for all subjects.
#' 
#' @details These are the details of the function
#' @return A list with some stuff
#' @export
#' @author Jonathan E. Gellar <jgellar1@@jhu.edu>
#' 

hf <- function(X, tind=NULL, basistype = c("s", "te", "t2"),
               additive=FALSE, concurrent=FALSE,
               integration = c("simpson", "trapezoidal", "riemann"),
               limits=NULL, splinepars=NULL) {
  tt.func <- create.tt.hf(X, tind=tind, basistype = basistype,
                          additive=additive, concurrent=concurrent,
                          integration = integration,
                          limits=limits, splinepars=splinepars)
  attr(X, "tt") <- tt.func
  X
}

hf <- function(X, ...) {
  tt.func <- create.tt.hf(X, ...)
  attr(X, "tt") <- tt.func
  X
}



hlf <- function(X, tind=NULL, basistype = c("s", "te", "t2"),
                integration = c("simpson", "trapezoidal", "riemann"),
                limits=NULL, splinepars=NULL) {
  tt.func <- create.tt.hlf()
  attr(X, "tt") <- tt.func
  X
}

haf <- function(X, sind=NULL, basistype = c("s", "te", "t2"), integration = c("simpson", "trapezoidal", "riemann"),
                limits=NULL, splinepars=NULL) {
  tt.func <- create.tt.haf()
  attr(X, "tt") <- tt.func
  X
  #   list(X=X, opts=list(sind=sind, basistype=basistype, integration=integration, limits=limits, splinepars=splinepars))
}
