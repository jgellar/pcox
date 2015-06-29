#' Create x-transform function for non-time-varying terms
#' 
#' The \code{xt} function processes data that does not involve time
#' transformations for inclusion in the \code{coxph} formula. Unlike a \code{tt}
#' function for time-transformed terms (see \code{\link{create.tt.func}}),
#' this function is not required by \code{coxph}. However, we use it for
#' consistency with how we handle \code{tt} terms, and also because it provides
#' a convenient place to store all information needed to process the data
#' and set up any required basis.
#' 
#' @param limits specifies the term as either a term involving scalar
#'   covariates or a functional covariate.
#'   Defaults to \code{NULL}, indicating scalar covariates. See Details.
#' @param linear if \code{FALSE}, covariates are included as nonlinear (smooth)
#'   effects, otherwise as a linear effect.
#' @param basistype specifies the basis constructor function (from the
#'   \code{mgcv} package) that is used to define a smooth term. Defaults to
#'   \code{\link[mgcv]{s}}, which is the only option allowed for smooths of
#'   only one argument. For smooths of multiple arguments,
#'   \code{\link[mgcv]{te}} or \code{\link[mgcv]{t2}} may (but don't have to)
#'   be used.
#' @param sind specifies the argument values for functional predictors. Can be
#'   entered as a vector of length \code{ncol(X)}, or a matrix of the same
#'   dimensions as \code{X} (for covariates measured on unequal grids).
#' @param integration method for numerical integration.
#' @param standardize if \code{TRUE}, the term is "standardized" by dividing
#'   by the width of integration.
# @param s.transform optional transformation function for the first variable
#   of the smooth. For functional/historical predictors, this is the variable
#   over which the integration takes place.
# @param t.transform optional transformation function for the time variable,
#   if it is one of the indices of the smooth
#' @param basisargs arguments for the function specified by \code{basistype},
#'   used to set up the basis and penalization
#' @param method method of optimization of the smoothing parameter, for
#'   penalized terms
#' @param eps tolerance level for the optimization of the smoothing parameter,
#'   for penalized terms
#' 
#' @return An \code{xt} function.
#' The \code{xt} function takes one argument, \code{x} (the data). According
#' to what is entered for the \code{limits} argument, the \code{xt} function
#' will return one of the following:
#' \enumerate{
#'   \item For terms that do not involve smooths, it simply returns the model
#'         matrix columns for this term
#'   \item For terms that do involve smooths, it calls \code{pcoxTerm} and
#'         returns the result.
#' }
#' 
#' Note that the \code{xt} function will contain all the information needed
#' to process the data in its environment. For this reason, \code{pcox} saves
#' the \code{xt} function as part of its return object, so that any methods
#' called on that object (e.g., \code{predict} or \code{coef}) know how to
#' process the data.
#' 
#' @seealso \code{\link[survival]{coxph}}, \code{\link{p}},
#'   \code{\link{pcoxTerm}}, \code{\link{create.tt.func}}
#' @keywords internal

create.xt.func <- function(limits, linear, basistype, sind, integration,
                           standardize, #s.transform,# t.transform,
                           basisargs, method, eps) {
  
  # Initialize: no smooth object created yet
  smooth <- NULL
  
  # Process limits argument: create appropriate processing function
  if (!is.null(limits))
    # limits is either "all", "full", or "baseline" - baseline function
    limits <- function(s,t) {s==s} # Will just spit out TRUE's
  smooth.flag <- !(is.null(limits) & linear) # No smooth if its a linear scalar
  
  xt.func <- function(x) {
    if (smooth.flag) {
      # Smooth involved: call pcoxTerm
      pcoxTerm(x, limits=limits, linear=linear, tv=FALSE,
               basistype=basistype, sind=sind,
               integration=integration, standardize=standardize,
               basisargs=basisargs, method=method, eps=eps,
               env=env, index=index, smooth=smooth)
    } else {
      # No smooth involved: just return the data (glorified identity function)
      x
    }
  }
  xt.func
}
