#' Define penalized or time-varying terms in a pcox formula
#' 
#' Function used to set up special terms in a pcox formula. These terms include
#' penalized terms and/or terms with time-varying coefficients or effects.
#' 
#' @param ... a list of variables that are the covariates used in the term, as
#'   well as possibly additional arguments that are passed onto the basis
#'   constructor defined by \code{basistype}.
#' @param limits specifies the term as either a term involving scalar
#'   covariates, a concurrent effect of a time-varying covariate, a baseline
#'   functional covariate, or a historical effect of a time-varying covariate.
#'   Defaults to \code{NULL}, indicating scalar covariates. See Details.
#' @param linear if \code{FALSE}, covariates are included as nonlinear (smooth)
#'   effects, otherwise as a linear effect.
#' @param tv if \code{TRUE}, makes the effect time-varying.
#' @param basistype specifies the basis constructor function (from the
#'   \code{mgcv} package) that is used to define a smooth term. Defaults to
#'   \code{\link[mgcv]{s}}, which is the only option allowed for smooths of
#'   only one argument. For smooths of multiple arguments (including t and s),
#'   \code{\link[mgcv]{te}} or \code{\link[mgcv]{t2}} may (but don't have to)
#'   be used.
#' @param sind specifies the argument values for functional and time-varying
#'   predictors. Can be entered as a vector of length \code{ncol(X)}, or a
#'   matrix of the same dimensions as \code{X} (for covariates measured on
#'   unequal grids).
#' @param integration method for numerical integration.
#' @param standardize if \code{TRUE}, the term is "standardized" by dividing
#'   by the width of integration.
#'   
#' @details The \code{limits} argument defines the type of term. Options include:
#'   \enumerate{
#'     \item Scalar terms: \code{NULL}
#'     \item Baseline functional predictors: a character string, any of
#'       \code{"baseline"}, \code{"all"}, or \code{"full"}, indicating to use
#'       the entire integration range at all times.
#'     \item Concurrent time-varying covariates: any of the character strings
#'       \code{"t"}, \code{"s==t"}, or \code{"s=t"} indicate to use the current
#'       value of the covariate to effect the hazard at time t. Alternatively,
#'       a non-negative number may be used to indicate to use a lagged version
#'       of the time-varying covariate, lagged by the entered amount of time.
#'       If the covariate value is not available at exactly that time, the last
#'       value carried forward is used. 
#'     \item Historical time-varying covariates: one of \code{"s<t"} or
#'       \code{"s<=t"} to indicate the range of integration up to (and possibly
#'       including) time t. Alternatively, could be a function of \code{s} and
#'       \code{t} (in that order), which returns \code{TRUE} if the covariate
#'       value at time \code{s} should impact the hazard at time \code{t},
#'       and \code{FALSE} othwerwise. Note that this flexibility could result
#'       in unpredictable results, so if you enter a function, use at your own
#'       risk!
#'   }
#'   
#'   More details to come....
#' 
#' @author Jonathan Gellar <jgellar1@@jhu.edu>
#' @export
#' @return The return object is a list with two components. The first is the
#'   raw data required for the term. The second is either a time-transform
#'   function or an x-transform function.
#'   
#'   Time-transform functions (element $tt) are returned if the term involves
#'   any time-varying components (i.e., concurrent or historical time-varying
#'   covariate, or time-varying effect). It is a function of \eqn{x} and \eqn{t}
#'   that specifies how to set up the term within \code{coxph()}.
#'   
#'   x-transform functions (element $xt) are returned when the term does not
#'   involve any time-varying aspects. This is a function of \eqn{x} that
#'   specifies how to set up the term within \code{pcox()}. For simple scalar
#'   terms, this is just an identity function, but for smooth terms, it requires
#'   setting up a basis matrix.
#' 
#' @seealso \code{\link{bf}}, \code{\link{cf}}, and \code{\link{hf}}, which
#'   are user-friendly wrappers for \code{p} that provide default arguments for
#'   baseline functional, concurrent, and historical terms, respectively.
#'   Also, \code{\link[mgcv]{s}}, \code{\link[mgcv]{te}}, and
#'   \code{\link[mgcv]{t2}} for options available for each \code{basistype},
#'   as well as \code{mgcv}'s \code{\link[mgcv]{smooth.terms}} for details of 
#'   \code{mgcv} syntax and available spline bases and penalties.
#' 

p <- function(..., limits=NULL, linear = TRUE, tv = FALSE,
              basistype = c("s", "te", "t2"), sind=NULL,
              integration = c("riemann", "trapezoidal", "simpson"),
              standardize = FALSE
              #domain = c("s", "s-t", "s/t"),
              ) {
  basistype <- match.arg(basistype)
  integration <- match.arg(integration)
  #domain <- match.arg(domain)
  
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
  
  # Determine if the term is time-varying or time-static, based on limits/tv
  tt <- if (is.null(limits)) {
    tv
  } else if (is.character(limits)) {
    if (tolower(limits) %in% c("all", "full", "baseline"))
      tv
    else TRUE
  } else if (is.function(limits)) {
    length(formals(limits)) > 1
  } else if (is.numeric(limits)) {
    TRUE
  } else stop("Unrecognized entry for limits argument")
  
  if (!tt) {
    # No time-varying aspect to the term: create a xt function
    xt <- create.xt.func(limits, linear, basistype, sind, integration,
                         standardize, basisargs, method, eps)
    # Return
    list(x=data, xt=xt)
  } else {
    # Time-varying aspect to the term: create a tt function
    
    # data must be passed as a vector or matrix (can't be a data frame)
    # The map will be used to recreate the data frame within the tt function
    map <- vector("list", length=length(data))
    cnt <- 0
    for (i in 1:length(data)) {
      nc <- ifelse(is.null(ncol(data[[i]])), 1, ncol(data[[i]]))
      map[[i]] <- (cnt+1):(cnt+nc)
      cnt <- cnt+nc
    }
    names(map) <- names(data)
    x <- as.matrix(data)
    
    # Create the tt function
    tt <- create.tt.func(limits, linear, tv, basistype, sind, integration,
                         standardize, basisargs, method, eps, map)
    # Return
    list(x=x, tt=tt)    
  }
}
