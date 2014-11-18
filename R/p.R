#' Define penalized or time-varying terms in a pcox formula
#' 
#' Function used to set up special terms in a pcox formula. These terms include
#' penalized terms and/or terms with time-varying coefficients or effects.
#' 
#' @param ... a list of variables that are the covariates used in the term, as well
#'   as possibly additional arguments that are passed onto the basis constructor
#'   defined by \code{basistype}
#' @param limits specifies the term as either a term involving scalar covariates,
#'   a concurrent effect of a time-varying covariate, a baseline functional covariate,
#'   or a historical effect of a time-varying covariate. Defaults to \code{NULL},
#'   indicating scalar covariates. See Details.
#' @param linear if \code{FALSE}, covariates are included as nonlinear (smooth)
#'   effects, otherwise as a linear effect
#' @param tv if \code{TRUE}, makes the effect time-varying
#' @param basistype specifies the basis constructor function (from the
#'   \code{mgcv} package) that is used to define a smooth term. Defaults to
#'   \code{\link[mgcv]{s}, which is the only option allowed for smooths of only
#'   one argument. For smooths of multiple arguments (including t and s), 
#'   \code{\link[mgcv]{te} or \code{\link[mgcv]{t2} may (but don't have to)
#'   be used.
#' @param sind specifies the time indices for functional and time-varying
#'   predictors. Can be entered as a vector of length \code{ncol(X)}, or a
#'   matrix of the same dimensions as \code{X} (for covariates measured on
#'   unequal grids).
#' 
#' @details These are the details... lots to go in here
#' 
#' 
#' @author Jonathan Gellar <jgellar1@@jhu.edu>
#' @return If the term involves time in any way, i.e., if it is a time-varying
#'   covariate (concurrent or historical) or time-varying effect, the return
#'   object is a list with two arguemnts: the raw data required for the term,
#'   and a function of \eqn{x} and \eqn{t} that specifies how to set up the
#'   term within \code{coxph()}. If the term is a smooth term that doesn't
#'   involve time, the return object is a list with two objects: a
#'   \code{coxph.penalty} object, and a list of the \code{smooth} objects
#'   from \code{mgcv::smoothCon()} that contain the basis information. If
#'   the term involves neither time nor a penalty, the data is simply
#'   returned.
#' @seealso \code{mgcv}'s \code{\link[mgcv]{smooth.terms}} for details of 
#'   \code{mgcv} syntax and available spline bases and penalties; the related 
#'   \code{\link[refund]{pffr}} and \code{\link[refund]{fgam}} from 
#'   \code{refund}.
#' 

p <- function(..., limits=NULL, linear = TRUE, tv = FALSE,
              basistype = c("s", "te", "t2"), sind=NULL,
              dbug=FALSE) {
  # Separate data from basis options
  frmls <- formals(match.fun(basistype))
  dots  <- list(...)
  args  <- names(dots) %in% names(frmls)
  basisargs <- dots[args]
  data      <- as.data.frame(lapply(dots[!args], I))
  names(data) <- as.list(substitute(list(...)))[-1][!args]
  
  if ((is.null(limits) | tolower(limits) %in% c("all", "full")) & !tv) {
    # No tt function required
    if (is.null(limits) & linear) {
      # No smooth required - just return data
      data
    } else {
      # Make coxph.penalty term via pcoxTerm
      if (!is.null(limits)) limits <- function(s,t) TRUE
      pcoxTerm(data, limits=limits, linear=linear, tv=tv,
               basistype=basistype, sind=sind, basisargs)
    }
  } else {
    # tt function required: create function
    tt <- create.tt.p(limits, linear, tv, basistype, sind, basisargs)
    list(data=data, tt=tt)
  }
}