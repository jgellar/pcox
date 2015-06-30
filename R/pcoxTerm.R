#' Create a smooth and coxph.penalty object for use in a pcox formula
#' 
#' This function, called by either a \code{tt} or \code{xt} function, will
#' create a new \code{mgcv}-style smooth object, save it, and convert it
#' into a \code{coxph.penalty} term (by calling \code{\link{pterm}}), as long
#' as the smooth object is not already created. If the smooth object has
#' already been created (and is supplied as the \code{smooth} argument),
#' then \code{pcoxTerm} returns the prediction matrix (by calling
#' \code{\link[mgcv]{PredictMat}}).
#' 
#' @param data data frame needed to set up the smooth term or prediction matrix
#' @param limits either a function to define integration limits for
#'   functional/historical predictors, or \code{NULL} indicating a smooth of
#'   scalar variables.
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
#' @param basisargs arguments for the function specified by \code{basistype},
#'   used to set up the basis and penalization
#' @param method method of optimization of the smoothing parameter, for
#'   penalized terms
#' @param eps tolerance level for the optimization of the smoothing parameter,
#'   for penalized terms#' 
#' @param env environment that contains the list of \code{smooth} objects
#'   within \code{pcox()}. This allows \code{pcoxTerm} to save a newly generated
#'   smooth somewhere that \code{pcox} has access to. This cannot be done by
#'   returning the \code{smooth} object for cases when \code{pcoxTerm} is called
#'   within a \code{tt} function by \code{coxph}.
#' @param index the index of \code{env$smooth} where a newly generated smooth
#'   should be stored
#' @param smooth an optional supplied \code{mgcv}-style smooth object. If it
#'   is present, then that object is used to generate a prediction matrix for
#'   the new \code{data}. If \code{NULL}, a new smooth object is created.
#' @param t time points for time-varying terms
#' 
#' @return One of the following:
#' \enumerate{
#'   \item If \code{smooth} is \code{NULL}, a \code{coxph.penalty} term
#'         resulting from a call to \code{\link{pterm}}
#'   \item If \code{smooth} is an \code{mgcv}-style smooth object, the
#'         prediction matrix resulting from a call to
#'         \code{\link[mgcv]{PredictMat}}.
#' }
#' @keywords internal
#' 

pcoxTerm <- function(data, limits, linear, tv, basistype, sind, integration,
                     standardize, basisargs, method, eps, env, index,
                     smooth, t=NULL) {
  
  data[[1]][is.na(data[[1]])] <- 0
  evaldat <- data
  newcall <- list(as.symbol(basistype))
  varnames <- lapply(names(data), as.symbol)
  n.var <- length(varnames)
  
  if (is.null(limits)) {
    # Smooth of scalar variable(s)
    if (tv) {
      evaldat$t <- t
      newcall <- c(newcall, quote(t))
    }
    if (!linear)
      newcall <- c(newcall, varnames)
    else if (n.var==1)
      newcall <- c(newcall, by=varnames)
    else {
      warning("Multiple smooth terms, using last one as by variable")
      newcall <- c(newcall, varnames[-n.var],
                   by=varnames[n.var])
    }
  } else {
    # Functional Term
    n <- nrow(data[[1]])
    J <- ncol(data[[1]])
    smat <- if (is.matrix(sind)) sind
    else matrix(sind, nrow=n, ncol=J, byrow=TRUE)
    
    if (!is.null(t)) {
      tmat <- matrix(t, nrow=n, ncol=J)
      mask <- t(outer(smat[1,], tmat[,1], limits))
    } else {
      # t will be NULL if it's not a tt term: assume full range
      mask <- matrix(TRUE, nrow=n, ncol=J)
      tmat <- NULL
    }
    L <- getL3(smat, integration, mask)
    if (standardize)
      L <- L / matrix(apply(smat*mask, 1, function(x) diff(range(x))),
                      nrow=nrow(L), ncol=ncol(L))
    
    # (Temporarily) replace masked out coordinates with NA
    smat[!mask] <- NA
    if (!is.null(t)) tmat[!mask] <- NA
    
    # Replace mased coordinates of smat and tmat
    smat[!mask] <- smat[mask][1]
    if (!is.null(t)) tmat[!mask] <- tmat[mask][1]
    
    smat.name <- paste0(varnames[1], ".smat")
    evaldat[[smat.name]] <- smat
    newcall <- c(newcall, as.symbol(substitute(smat.name)))
    if (tv) {
      # Allow smooth to vary over t
      tmat.name <- paste0(varnames[1], ".tmat")
      evaldat[[tmat.name]] <- tmat
      newcall <- c(newcall, as.symbol(substitute(tmat.name)))
    }
    
    if (!linear) {
      # Nonlinear functional term: X in smooth, by=L
      # possible transformations of X?
      #  Quantile
      #  Log
      #  Probit
      L.name <- paste0(varnames[1], ".L")
      evaldat[[L.name]] <- L
      newcall <- c(newcall, varnames, by=as.symbol(substitute(L.name)))
    } else if (length(data)==1) {
      # Linear functional term: by=LX
      LX <- data[[1]]*L
      LX.name <- paste0(varnames[1], ".LX")
      evaldat[[LX.name]] <- LX
      newcall <- c(newcall, by=as.symbol(substitute(LX.name)))
    } else {
      # Linear smooth with multiple functions...
      stop("Not sure how to handle smooths with multiple linear functions yet")
    }
  }
  
  # Return is based on whether or not a smooth object has been supplied.
  # If not, we create the new smooth and the coxph.penalty object, and
  # return them. If smooth is supplied, return prediction matrix
  if (is.null(smooth)) {
    # Createsmooth object
    newcall <- c(newcall, basisargs)
    smooth <- mgcv::smoothCon(eval(as.call(newcall)), data=evaldat, knots=NULL,
                        absorb.cons=TRUE)
    if (length(smooth)>1) {
      # Can we turn these into a single smooth object (wider basis matrix, 
      #   block diagonal penalty matrix)?
      stop("We don't yet support terms with multiple smooth objects.")
    }
    
    # Assign smooth object and smoothdata back to pcox
    env$smooth[[index]] <- smooth
    env$smoothdata[[index]] <- evaldat
    
    # Create coxph.penalty term based on smooth
    pterm(smooth[[1]], method=method, eps=eps)
  } else {
    # Return prediction matrix
    mgcv::PredictMat(smooth[[1]], data=evaldat)
  }
}