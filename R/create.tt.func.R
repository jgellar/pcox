#' Create \code{tt} function for a time-varying term
#' 
#' A time-transform (\code{tt}) function is used by \code{survival::coxph} to
#' allow variables to vary dynamically in time. This allows both time-varying
#' covariates and time-varying effects to be included in a \code{coxph} (and
#' thus \code{pcox}) formula. This function sets up the \code{tt} function
#' according to the supplied arguments, and returns that function for later
#' use by \code{coxph}.
#' 
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
#' @param sind specifies the time indices for functional and time-varying
#'   predictors. Can be entered as a vector of length \code{ncol(X)}, or a
#'   matrix of the same dimensions as \code{X} (for covariates measured on
#'   unequal grids).
#' @param integration method for numerical integration.
#' @param standardize if \code{TRUE}, the term is "standardized" by dividing
#'   by the width of integration.
#' @param s.transform optional transformation function for the first variable
#'   of the smooth. For functional/historical predictors, this is the variable
#'   over which the integration takes place.
#' @param t.transform optional transformation function for the time variable,
#'   if it is one of the indices of the smooth
#' @param basisargs arguments for the function specified by \code{basistype},
#'   used to set up the basis and penalization
#' @param method method of optimization of the smoothing parameter, for
#'   penalized terms
#' @param eps tolerance level for the optimization of the smoothing parameter,
#'   for penalized terms
#' @param map variable map used to turn the \code{x} argument of the \code{tt}
#'   function from a matrix back into a data frame. Necessary because \code{x}
#'   must be a matrix, but we need to know which columns correspond to which
#'   variables.
#' 
#' @return a \code{tt} function.
#' The \code{tt} function takes two arguments, \code{x} (the data) and \code{t}
#' (the time vector). For terms involving time-varying effects or covariates,
#' these arguments will already be expanded (by \code{coxph}) to "long" format
#' at the time the \code{tt} function is called. According to what is entered
#' for the \code{limits} argument, the \code{tt} function will do one of the
#' following:
#' \enumerate{
#'   \item For terms that do not involve smooths, it simply returns the model
#'         matrix columns for this term
#'   \item For terms that do involve smooths, it calls \code{pcoxTerm} and
#'         returns the result.
#' }
#' 
#' Note that the \code{tt} function will contain all the information needed
#' to process the data in its environment. For this reason, \code{pcox} saves
#' the \code{tt} function as part of its return object, so that any methods
#' called on that object (e.g., \code{predict} or \code{coef}) know how to
#' process the data.
#' 
#' @seealso \code{\link[survival]{coxph}}, \code{\link{p}},
#'   \code{\link{pcoxTerm}}, \code{\link{create.xt.func}}
#' @keywords internal

create.tt.func <- function(limits, linear, tv, basistype, sind, integration,
                           standardize, s.transform, t.transform,
                           basisargs, method, eps, map) {
  
  # Initialize: no smooth object created yet
  smooth = s0 = t0 <- NULL
  
  # Process limits argument: create appropriate processing function
  conc.fcn <- NULL
  if (is.character(limits)) {
    if (limits %in% c("s=t", "s==t", "t")) {
      #conc.fcn <- function(s,t) s==t
      conc.fcn <- function(s,t) abs(s-t)
      limits <- NULL
    } else if (tolower(limits) %in% c("all", "full")) {
      limits <- function(s,t) {s==s}
    } else if (limits=="s<t") {
      limits <- function(s,t) s<t
    } else if (limits=="s<=t") {
      limits <- function(s,t) s<=t
    } else {
      stop("Unrecognized limits option!")
    }
  } else if (is.numeric(limits)) {
    if (limits<0)
      stop("limits must be non-negative if using it to specify lagged 
           time-varying covariates")
    lag <- limits
    conc.fcn <- function(t, s) {
      targ <- t-lag
      ifelse(targ<min(s), NA, which.min(abs(s-targ)))
    }
    #conc.fcn <- function(s,t) abs(s - (t-lag))
    limits <- NULL # Treat the lagged covariate as a scalar
  } else if (!is.function(limits) & !is.null(limits)) {
    stop("Unrecognized input for limits argument")
  }
  
  # Process optional t transformation
  t.transform <- if (is.null(t.transform)) {
    if (!is.null(s.transform) & basistype=="s") {
      if (s.transform %in% c("s/t", "range")) # Defaults to "tmaxmin"
        function(t) (t-min(t0, na.rm=T))/(max(t0, na.rm=T)-min(t0, na.rm=T))
      else NULL            # Defaults to no transform
    } else NULL            # Defaults to no transform
  } else if (is.character(t.transform)) {
    if (t.transform=="t") NULL
    else if (t.transform=="tmax") function(t) t/max(t0, na.rm=T)
    else if (t.transform=="tmaxmin") function(t)
      (t-min(t, na.rm=T))/(max(t0, na.rm=T)-min(t0, na.rm=T))
    else if (t.transform=="ecdf") function(t) {
      y <- ecdf(t0)(t)
      if (is.matrix(t)) matrix(y, nrow=nrow(t), ncol=ncol(t))
      else y
    }
    else if (t.transform=="ecdf2") function(t) {
      stop("not implemented yet")
      y <- ecdf(t0)(t)
      y <- (y-min(y))/(max(y)-min(y))
      if (is.matrix(t)) matrix(y, nrow=nrow(t), ncol=ncol(t))
      else y
    }
    else stop("Unrecognized t transformation")
  } else if (!is.function(t.transform))
    stop("Unrecognized t tranformation: must be a function or a
         recognized transformation string")
  else if (length(formals(t.transform))>2 )
    stop("t.transform can only have 1 or 2 arguments")
  
  # Process optional s transformation
  s.transform <- if (is.null(s.transform))
    # Defaults to no transform
    NULL
  else if (is.character(s.transform)) {
    if (s.transform=="s") NULL
    else if (s.transform=="s-t") function(s,t) s-t
    else if (s.transform=="s/t") function(s,t) {
      smin <- min(s0, na.rm=TRUE)
      ifelse(t==smin, 0.5, (s-smin)/(t-smin))} 
    else if (s.transform=="range") function(s,t) {
      smt <- s-t
      t(apply(smt, 1, function(x) {
        y <- x[!is.na(x)]
        x[!is.na(x)] <- if (length(y)==1) 0.5 else (y-min(y))/(max(y)-min(y))
        x
      }))
    }
    else stop("Unrecognized s transformation")
  } else if (!is.function(s.transform))
    stop("Unrecognized s tranformation: must be a function or a
         recognized transformation string")
  else if (length(formals(s.transform))>2 | length(formals(s.transform))<1)
    stop("s.transform must have 1 or 2 arguments")
  
  tt.func <- function(x, t, ...) {
    
    # Turn x matrix back into a data.frame according to map
    if (is.vector(x)) {
      data <- data.frame(x)
    } else {
      rownames(x) <- NULL
      data <- as.data.frame(lapply(map, function(map.i) {
        if (length(map.i>1))
          I(x[,map.i])
        else
          x[,map.i]
      }))
    }
    names(data) <- names(map)
    
    # Process concurrent time-varying covariates
    if (!is.null(conc.fcn)) {
      # Reduce data matrices to vectors
      cfidx <- sapply(data, is.matrix)
      data[cfidx] <- lapply(data[cfidx], function(x) {
        if (ncol(x)!=length(sind))
          stop("Mismatch between length of sind and number of columns 
               of data matrix")
        #idxs <- sapply(t, function(t.i) which.min(conc.fcn(sind, t.i)))
        #sapply(1:nrow(x), function(i) x[i,idxs[i]])
        idxs <- sapply(t, conc.fcn, s=sind)
        x[cbind(1:nrow(x), idxs)]
      })
      names(data)[cfidx] <- paste0(names(data)[cfidx], ".t")
    }
    
    if (is.null(limits) & !tv & linear) {
      # No penalized term required: just return data 
      #   (i.e., the most basic concurrent TVC - no extras)
      as.matrix(data)
    } else {
      # Create and return pcoxTerm
      pcoxTerm(data, limits, linear, tv, basistype, sind, integration,
               standardize, s.transform, t.transform,
               basisargs, method, eps, env, index,
               smooth, s0, t0, t)
    }
  }
  # Return
  tt.func
}
