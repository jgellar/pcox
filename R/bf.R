#' Define baseline functional terms in a pcox formula
#' 
#' Function used to define a baseline functional term in a pcox formula.
#' Simply calls \code{\link{p}} with default arguments for baseline functions.
#' 
#' @param ... a list of variables that are the covariates used in the term, as well
#'   as possibly additional arguments that are passed onto the basis constructor
#'   defined by \code{basistype}
#' @param linear if \code{FALSE}, covariates are included as nonlinear (smooth)
#'   effects, otherwise as a linear effect
#' @param tv if \code{TRUE}, makes the effect time-varying
#' @param basistype character string that specifies the basis constructor
#'   function (from the \code{mgcv} package) that is used to define a smooth
#'   term. For linear, non-time varying functional effects, must be \code{"s"}
#'   because the smooth will only be over one argument. For other types of
#'   functional effects, \code{"te"} or \code{"t2"} may (but don't have to)
#'   be used.
#' @param sind specifies the time indices for functional predictor. Can be
#'   entered as a vector of length \code{ncol(X)}, or a matrix of the same
#'   dimensions as \code{X} (for covariates measured on unequal grids).
#'   Defaults to an equally spaced grid from 0 to 1.
#' @param integration method for numerical integration over \code{sind}
#' 
#' @details Baseline functional effects estimate a weight function across
#'   the domain of the functional predictor, and integrate over this effect
#'   to get the total contribution of the predictor towards the (log) hazard
#'   function. These effects differ from historical functional effects
#'   because the domain of the function here is assumed to not be related to
#'   the study time \eqn{t}.
#'   
#'   The default baseline functional effect is linear and non-time-varying: 
#'   \eqn{\int_0^S X_i(s)\beta(s) ds}. A nonlinear effect will move the
#'   functional predictor inside the \eqn{\beta} function, i.e.,
#'   \eqn{\int_0^S \beta[ s, X_i(s)] ds}. A time-varying effect will add
#'   \eqn{t} to the arguments of the \eqn{\beta} function, i.e.
#'   \eqn{\int_0^S X_i(s)\beta(s, t) ds} or
#'   \eqn{\int_0^S \beta[ s,t,X_i(s)] ds}.
#' 
#' @author Jonathan Gellar <jgellar1@@jhu.edu>
#' @export
#' @return If the term is time-varying, the return object is a list with
#'   two arguemnts: the raw data required for the term, and a function of
#'   \eqn{x} and \eqn{t} that specifies how to set up the term within
#'   \code{coxph()}. Otherwise, the return object is a list with two objects:
#'   a \code{coxph.penalty} object, and a list of the \code{smooth} objects
#'   from \code{mgcv::smoothCon()} that contain the basis information.
#'   
#' @seealso \code{\link{p}}
#' 

bf <- function(..., linear = TRUE, tv = FALSE, basistype = c("s", "te", "t2"),
               sind=NULL, integration=c("riemann", "trapezoidal", "simpson")#,
               #vd=FALSE
               ) {
  
  # Check basistype against the number of arguemnts of the smooth?
  basistype <- match.arg(basistype)
  if (linear & !tv & basistype!="s") {
    warning("Smooths over a single argument must use basistype \"s\"")
    basistype <- "s"
  }
  if (is.null(sind)) {
    dots <- list(...)
    dat1 <- dots[names(dots)==""][[1]]
    sind <- seq(0, 1, length=ncol(dat1))
  }
  
  p(..., limits="all", linear=linear, tv=tv, basistype=basistype, sind=sind,
    integration=integration, standardize=FALSE)
}
