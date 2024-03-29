#' Define historical (functional) time-varying covariate term in a pcox formula
#' 
#' Function used to define a historical functional term in a pcox formula.
#' Simply calls \code{\link{p}} with default arguments for historical terms.
#' 
#' @param ... a list of variables that are the covariates used in the term, as well
#'   as possibly additional arguments that are passed onto the basis constructor
#'   defined by \code{basistype}. Historical time-varying covariates should be
#'   included as an \eqn{N x J} matrix, where \eqn{N} is the number of subjects
#'   and \eqn{J} is the number of time points. The right-most columns of most rows
#'   will have \code{NA} values if they are after the event/censoring time.
#' @param limits specifies the range of integration for the historical effect.
#'   May be a function of \eqn{(s,t)} that returns \code{TRUE}/\code{FALSE} based
#'   on whether or not on observation at that combination of \eqn{(s,t)} should or
#'   should not be included. Alternatively, the character strings \code{"s<=t"}
#'   or \code{"s<t"} may be entered.
#' @param linear if \code{FALSE}, covariates are included as nonlinear (smooth)
#'   effects, otherwise as a linear effect
#' @param tv if \code{TRUE}, makes the effect time-varying
#' @param basistype character string that specifies the basis constructor
#'   function (from the \code{mgcv} package) that is used to define a smooth
#'   term. For linear, non-time varying functional effects, must be \code{"s"}
#'   because the smooth will only be over one argument. For other types of
#'   functional effects, \code{"te"} or \code{"t2"} may (but don't have to)
#'   be used.
#' @param sind specifies the time indices for the time-varying covariate. May
#'   be entered as a vector of length \code{ncol(X)}, or a matrix of the same
#'   dimensions as \code{X} (for covariates measured on unequal grids).
#' @param integration method for numerical integration over \code{sind}
#' @param standardize standardize term by dividing by the integration width?
#' @param transform character string indicating an optional basis transformation;
#'    see Details for options.
#' 
#' 
#' @details Historical functional effects involve time-varying covariates.
#'   They differ from concurrent effects in that for a historical effect,
#'   the entire (or partial) history of the time-varying covariate (up to
#'   time \eqn{t}) is allowed to effect the (log) hazard function at time
#'   \eqn{t}. This is accomplished by estimating a weight function across
#'   the time domain, and integrating this effect to get the total contribution
#'   of the time-varying covariate history towards the (log) hazard function
#'   at time \eqn{t}. These effects differ from baseline functional effects
#'   because in the latter, the domain of the function is independent of
#'   \eqn{t}.
#'   
#'   The default integration range is over the entire history of the covariate,
#'   i.e., from 0 to \eqn{t}. Unlike other types of \code{pcox} terms, the
#'   default is to assume a time-varying historical effect:
#'   \eqn{\int_0^t X_i(s)\beta(s,t) ds}. This may be overridden by specifying
#'   \code{tv=FALSE}. The term may also be made nonlinear, i.e.
#'   \eqn{\int_0^t \beta [s,t,X_i(s)] ds} or
#'   \eqn{\int_0^t \beta [s,X_i(s)] ds}.
#' 
#' @author Jonathan Gellar <jgellar1@@jhu.edu>
#' @export
#' @return The result of a call to \code{p()}, which will be a list with
#'   the raw data required for the term, and a function of \eqn{x} and \eqn{t}
#'   that specifies how to set up the term within \code{coxph()}.
#' @seealso \code{\link{p}}
#' 

hf <- function(..., limits = "s<=t", linear = TRUE, tv = TRUE,
               basistype = c("s", "te", "t2"), sind=NULL,
               integration=c("riemann", "trapezoidal", "simpson"),
               standardize = TRUE, transform = NULL) {
  basistype <- match.arg(basistype)
  integration <- match.arg(integration)
  #domain <- match.arg(domain)
  
  # Do some checks?
  if (is.null(sind)) {
    dots <- list(...)
    mat1 <- which(sapply(dots, is.matrix))[[1]]
    sind <- 1:ncol(dots[[mat1]])
  }
  
  # Process transformation
  if (!is.null(transform)) {
    # Set up new call to p, with bs and xt updated
    localP <- function(..., bs=NULL, xt=NULL, mp=NULL) {
      # Set up xt info for "dt" basis
      newxt <- switch(transform,
                      lagged = list(tf=list("s-t"), bs=bs, xt=xt),
                      standardized = list(tf=list("s/t", "linear01"), bs=bs, xt=xt),
                      noInteraction = list(tf="s/t", bs="pi",
                                           xt=list(g="none", bs=bs, xt=xt)),
                      linear = list(tf=list("s/t", "linear01"), bs="pi",
                                    xt=list(g="linear", bs=bs, xt=xt, mp=mp)),
                      quadratic = list(tf=list("s/t", "linear01"), bs="pi",
                                       xt=list(g="quadratic", bs=bs, xt=xt, mp=mp))
      )
      newxt$basistype <- basistype
      newxt <- rmNullObs(newxt)
      
      p(..., limits=limits, linear=linear, tv=tv, basistype="s", sind=sind,
        integration=integration, standardize=standardize,
        bs="dt", xt=newxt)
    }
    localP(...)
  } else {
    # Call p directly
    p(..., limits=limits, linear=linear, tv=tv, basistype=basistype, sind=sind,
      integration=integration, standardize=standardize)
  }
}

is.NullOb <- function(x) is.null(x) | all(sapply(x, is.null))

## Recursively step down into list, removing all such objects 
rmNullObs <- function(x) {
  x <- Filter(Negate(is.NullOb), x)
  lapply(x, function(x) if (is.list(x)) rmNullObs(x) else x)
}

