#' Define concurrent time-varying covariate term in a pcox formula
#' 
#' Function used to define a concurrent time-varying covariate term in a
#' pcox formula. Simply calls \code{\link{p}} with default arguments for
#' concurrent terms.
#' 
#' @param ... a list of variables that are the covariates used in the term, as well
#'   as possibly additional arguments that are passed onto the basis constructor
#'   defined by \code{basistype}. Concurrent time-varying covariates should be
#'   included as an \eqn{N x J} matrix, where \eqn{N} is the number of subjects
#'   and \eqn{J} is the number of time points. The right-most columns of most rows
#'   will have \code{NA} values if they are after the event/censoring time.
#' @param lag optional time lag for the concurrent effect. Term will be processed
#'   as \eqn{X_i(t-lag)}. Defaults to 0 (no lag).
#' @param linear if \code{FALSE}, covariates are included as nonlinear (smooth)
#'   effects, otherwise as a linear effect
#' @param tv if \code{TRUE}, makes the effect time-varying
#' @param basistype character string that specifies the basis constructor
#'   function (from the \code{mgcv} package) that is used to define a smooth
#'   term. For concurrent terms this is only relevant if the term is
#'   time-varying and nonlinear, because smooths of one variable must use
#'   \code{s}.
#' @param sind specifies the time indices for the time-varying covariate. May
#'   be entered as a vector of length \code{ncol(X)}, or a matrix of the same
#'   dimensions as \code{X} (for covariates measured on unequal grids).
#' 
#' @details A concurrent functional term in a \code{pcox} formula is a term
#'   that involves \eqn{X_i(t)}, or alternatively \eqn{X_i(t-lag)}. This term
#'   may have a linear or nonlinear effect of \eqn{X_i(t)} and may be time-fixed
#'   or time-varying (e.g., \eqn{\beta(t)X_i(t)}).
#' 
#' @author Jonathan Gellar <jgellar1@@jhu.edu>
#' @export
#' @return The result of a call to \code{p()}, which will be a list with
#'   the raw data required for the term, and a function of \eqn{x} and \eqn{t}
#'   that specifies how to set up the term within \code{coxph()}.
#'   
#' @seealso \code{\link{p}}
#' 

cf <- function(..., lag=0, linear = TRUE, tv = FALSE,
               basistype = c("s", "te", "t2"),
               sind=NULL, dbug=FALSE) {
  
  if (is.null(sind)) {
    dots <- list(...)
    sind <- 1:ncol(dots[names(dots)==""][[1]])
  }
  
  p(..., limits=(-lag), linear=linear, tv=tv, basistype=basistype, sind=sind,
    standardize=FALSE, s.transform=NULL, t.transform=NULL, dbug=dbug)
}

