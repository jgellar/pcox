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
#' @param integration method for numerical integration
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
#' @seealso \code{\link{bf}}, \code{\link{cf}}, and \code{\link{hf}}, which
#'   are wrappers for \code{p} that provide correct default arguments for
#'   baseline functional, concurrent, and historical terms, respectively.
#'   Also, \code{\link[mgcv]{s}}, \code{\link[mgcv]{te}}, and
#'   \code{\link[mgcv]{t2}} for options available for each \code{basistype},
#'   as well as \code{mgcv}'s \code{\link[mgcv]{smooth.terms}} for details of 
#'   \code{mgcv} syntax and available spline bases and penalties.
#' 

p <- function(..., limits=NULL, linear = TRUE, tv = FALSE,
              basistype = c("s", "te", "t2"), sind=NULL,
              integration=c("riemann", "trapezoidal", "simpson"),
              divide.by.t=FALSE, domain=c("s", "s-t", "u"), dbug=FALSE) {
  basistype <- match.arg(basistype)
  integration <- match.arg(integration)
  domain <- match.arg(domain)
  
  # Extract basis options
  frmls <- formals(match.fun(basistype))
  dots <- list(...)
  args <- if(is.null(names(dots)))
    rep(FALSE, length(dots))
  else  names(dots) %in% names(frmls)
  basisargs <- dots[args]
  
  # Extract method and eps
  method <- dots$method
  eps <- dots$eps
  
  # Set up data
  vars <- !args & !(names(dots) %in% c("method", "eps"))
  data      <- as.data.frame(lapply(dots[vars], I))
  names(data) <- as.list(substitute(list(...)))[-1][vars]
  
  if (any(is.null(limits), tolower(limits) %in% c("all", "full")) & !tv) {
    #if ((is.null(limits) | tolower(limits) %in% c("all", "full")) & !tv) {
    # No tt function required
    if (is.null(limits) & linear) {
      # No smooth required - just return data
      data
    } else {
      # Make coxph.penalty term via pcoxTerm
      if (!is.null(limits)) limits <- function(s,t) s==s
      pcoxTerm(data, limits=limits, linear=linear, tv=tv,
               basistype=basistype, sind=sind,
               integration=integration, divide.by.t=divide.by.t,
               domain=domain, basisargs=basisargs,
               method=method, eps=eps)
    }
  } else {
    # tt function required
    
    # First, create the data map and convert data to a matrix
    map <- vector("list", length=length(data))
    cnt <- 0
    for (i in 1:length(data)) {
      nc <- ifelse(is.null(ncol(data[[i]])), 1, ncol(data[[i]]))
      map[[i]] <- (cnt+1):(cnt+nc)
      cnt <- cnt+nc
    }
    names(map) <- names(data)
    x <- as.matrix(data)
    #attr(x, "map") <- map
    
    # Create the tt function
    tt <- create.tt.p(limits, linear, tv, basistype, sind, basisargs,
                      method=method, eps=eps, map=map)
    if (dbug) debug(tt)
    
    # data must be a vector or matrix (can't be a data.frame)
    list(x=x, tt=tt)
  }
}
