#' Extract coefficient functions from a fitted pcox-object
#' 
#' This function is used to extract coefficients from a fitted `pcox` model, in
#' particular smooth functions resulting from including smooth effects
#' of scalar covariates, time-varying effects, functional predictors, or
#' historical effects. The function calls \code{mgcv::plot.gam} to generate
#' the smooth functions.
#' 
#' @param object return object from \code{\link{pcox}}
#' @param select integer vector indicating the index of the desired smooth term
#'   in \code{object$pcox$smooth}. If \code{NULL}, returns all smooth terms.
#' @param se if \code{TRUE}, returns pointwise standard error estimates. If
#'   numeric, the standard error will be multiplied by that number (e.g., 2
#'   to make it easier to plot confidence intervals).
#' @param use.var2 if \code{TRUE}, uses \code{object$var2} as the covariance
#'   matrix for standard error estimation; otherwise uses \code{object$var}. 
#' @param exclude if \code{TRUE}, excludes reporting of the estimate at
#'   coordinates that are "too far" from data used to fit the model, as
#'   determined by \code{mgcv::plot.mgcv.smooth}. Excluding these values sets
#'   them to \code{NA} in the return object.
#' @param ... further arguments passed on to \pkg{mgcv}'s \code{plot.gam}. Common
#'   arguments include \code{n} and \code{n2} to set the number of coordinates
#'   to estimate for 1-D and 2-D coefficient functions, and \code{seWithMean}
#'   if the standard error should include uncertainty about the overall mean.
#'   See \code{\link[mgcv]{plot.gam}}.
#' 
#' @return a list of data.frames containing the evaluation points, 
#'    coefficient function values and optionally their se's for each term in \code{select}.
#'    If only one term is selected, the one data frame is unlisted.
#' 
#' @author Fabian Scheipl and Jonathan Gellar, adapted from
#'   \code{refund::coefficients.pfr}
# @examples 
# #TODO: see ?pfr 
#' @export

coefficients.pfr <- function(object, use.var2=FALSE, exclude=FALSE, ...) {
  
  # Create fakegam object
  fakegam <- list(coefficients=object$coefficients, cmX=object$means,
                  smooth=object$pcox$smooth, model=object$pcox$smoothdata[[1]])
  
  # Function to strip some items out of dots (...)
  localplot <- function(..., residuals=NULL, all.terms=NULL, too.far=0.1) {
    if (all(!is.null(residuals), residuals))
      warning("residuals option not allowed for coefficients.pfr; removing")
    if (all(!is.null(all.terms), all.terms))
      warning("all.terms option not allowed for coefficients.pfr; removing")
    if (!exclude) too.far <- 0
    if (any(is.null(list(...)$se), list(...)$se))
      fakegam$Vp <- if (use.var2) object$var2 else object$var
    
    mgcv::plot.gam(fakegam, too.far=too.far, ...)
  }
  
  plotdata <- localplot(...)
  
  #plotdata <- mgcv::plot.gam(object, residuals=FALSE, select=which, se=se,
  #                           too.far=too.far, ...)
  
  coef <- lapply(1:length(plotdata), function(i) {
    pd <- plotdata[[i]]
    
    # Create coef.i with coordinates
    is.re <- "random.effect" %in% class(object$smooth[[i]])
    if(is.re) pd$x <- levels(object$model[[object$smooth[[i]]$term]])
    
    coef.i <- if(object$smooth[[i]]$dim == 1) {
      setNames(data.frame(pd$x),
               ifelse(is.re, object$smooth[[i]]$term,
                      gsub("tmat", "argvals", pd$xlab)))
    } else {
      grid <- expand.grid(x=pd$x, y=pd$y)
      setNames(data.frame(grid$x, grid$y), 
               c(gsub("\\.omat", "", pd$xlab),
                 gsub("tmat", "argvals", pd$ylab)))
    }
    
    # Add estimate & SE
    coef.i$value <- pd$fit
    if (se)
      coef.i$se  <- pd$se
    
    coef.i
  })
  
  #special extra wish from Jon:
  if(length(coef) == 1) {
    coef[[1]]
  } else {
    names(coef) <- sapply(object$smooth, function(x) x$label)
    coef
  }
}
#' @rdname coefficients.pfr
coef.pfr <- coefficients.pfr


coefficients.pcox2 <- function(object, raw=FALSE, term=NULL, n=NULL, n2=NULL,
                               se=TRUE, seWithMean=TRUE, exclude=FALSE,
                               limit=TRUE, ...) {
  
  # Create fakegam object
  fakegam <- list(coefficients=object$coefficients, Vp=object$var,
                  smooth=object$pcox$smooth, model=object$pcox$smoothdata[[1]])
  tmp <- mgcv::plot.gam(fakegam)
  
  
  x$residuals - "partial.resids"
  x$pterms - "order"
  x$edf - same length as coefs, estimated df per term
  x$model - the data ******
  x$cmX - for seWithMean
  
}