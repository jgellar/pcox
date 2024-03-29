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
#' @param n number of estimation points for each 1-dimensional coefficient
#'   function. If estimating more than one coefficient, can be entered
#'   as a vector of length equal to the number of coefficients. Default
#'   is the number of unique observation points.
#' @param n2 Square root of the number of points used to grid estimates of
#'   2-dimensinoal functions for contouring.
#'   If estimating more than one coefficient, can be entered
#'   as a vector of length equal to the number of coefficients.
#' @param exclude if \code{TRUE}, excludes reporting of the estimate at
#'   coordinates that are "too far" from data used to fit the model, as
#'   determined by \code{mgcv::plot.mgcv.smooth}. Excluding these values sets
#'   them to \code{NA} in the return object.
#' @param limit if \code{TRUE}, checks if a \code{limits} function was
#'   used to generate the term, and if so, applies the function to the
#'   output to only produce estimates over regions that were within "limits".
#' @param plotMe if \code{TRUE}, also plots the coefficient(s)
#' @param ... further arguments passed on to \pkg{mgcv}'s \code{plot.gam}; see
#'   Details.
#'   
#' @details
#'   The smooth coefficients are generated by \code{mgcv::plot.gam}. The
#'   following are common arguments of \code{plot.gam} that can be
#'   specified in the function call to \code{coefficients.pcox}:
#'   \describe{
#'     \item{select}{integer vector indicating the index of the desired smooth
#'     term in \code{object$pcox$smooth}. If \code{NULL} (the default), returns
#'     all smooth terms.}
#'     \item{se}{if \code{TRUE} (the default), returns pointwise standard error
#'     estimates. The covariance matrix to use (either \code{var} or
#'     \code{var2}) is specified by \code{use.var2}.}
#'     \item{n}{number of evaluation points for each 1-D plot. Only equally-spaced
#'     grids are possible.}
#'     \item{n2}{Square root of the number of points used to grid estimates of
#'     2-D functions.}
#'   }
#'   See \code{\link[mgcv]{plot.gam}} for a full list of possible arguments.
#'   Note that many of these have to do with plotting, and will thus not affect
#'   the resultant coefficient function.
#'   
#'   This function currently does not work to extract coefficient involving
#'   smooths over 3 or more variables. To extract these coefficients, you
#'   must do so manually, see \code{\link[mgcv]{PredictMat}}.
#' 
#' @return a list of data.frames containing the evaluation points, 
#'    coefficient function values and optionally their se's for each term in \code{select}.
#'    If only one term is selected, the one data frame is unlisted and the return
#'    object will be the single data frame.
#' 
#' @author Fabian Scheipl and Jonathan Gellar, adapted from
#'   \code{refund::coefficients.pfr}
# @examples 
# #TODO: see ?pfr 
#' @export

coefficients.pcox2 <- function(object, select=NULL, se=TRUE, use.var2=FALSE,
                               n=NULL, n2=NULL,
                               exclude=FALSE, limit=TRUE, plotMe=FALSE, ...) {
  
  # Create fakegam object
  fakegam <- list(coefficients=object$coefficients, cmX=object$means,
                  smooth=object$pcox$smooth, model=object$pcox$smoothdata[[1]])
  if (se)
    fakegam$Vp <- if (use.var2) object$var2 else object$var
  
  # Function to strip some items out of dots (...)
  localplot <- function(..., residuals=NULL, all.terms=NULL, too.far=0.1) {
    if (all(!is.null(residuals), residuals))
      warning("residuals option not allowed for coefficients.pfr; removing")
    if (all(!is.null(all.terms), all.terms))
      warning("all.terms option not allowed for coefficients.pfr; removing")
    if (!exclude) too.far <- 0
    
    if (!plotMe) {
      ## dump plots to file since can't pass type = "n" to plot.gam
      tfile <- tempfile()
      pdf(tfile)
      plotdata <- mgcv::plot.gam(fakegam, too.far=too.far, ...)
      dev.off()
      if (file.exists(tfile))
        file.remove(tfile)
    } else
      plotdata <- mgcv::plot.gam(fakegam, too.far=too.far, ...)
    
    plotdata
  }
  
  if (is.null(select)) select <- seq_along(object$pcox$smooth)
  if (!is.null(n)) {
    if (length(n) ==1) n  <- rep(n,  length(select)) else if (length(n) !=length(select))
      stop("Mismatch between length of n and number of terms")
  }
  if (!is.null(n2)) {
    if (length(n2)==1) n2 <- rep(n2, length(select)) else if (length(n2)!=length(select))
      stop("Mismatch between length of n2 and number of terms")
  }
  
  coef <- lapply(select, function(i) {
    # Set values of n and n2 (if not supplied)
    sdat.i <- object$pcox$smoothdata[[i]]
    smooth.i <- object$pcox$smooth[[i]]
    n.i  <- ifelse(is.null(n[i]),  ndefault(sdat.i[[smooth.i$term[[1]]]]), n[i])
    n2.i <- ifelse(is.null(n2[i]), n.i, n2[i])
    
    pd <- localplot(..., n=n.i, n2=n2.i, select=i, se=se)[[1]]
    
    # Create coef.i with coordinates
    is.re <- "random.effect" %in% class(fakegam$smooth[[i]])
    if(is.re) pd$x <- levels(fakegam$model[[fakegam$smooth[[i]]$term]])
    
    coef.i <- if(fakegam$smooth[[i]]$dim == 1) {
      setNames(data.frame(pd$x),
               ifelse(is.re, fakegam$smooth[[i]]$term, modify_st(pd$xlab)))
    } else {
      grid <- expand.grid(x=pd$x, y=pd$y)
      setNames(data.frame(grid$x, grid$y),
               c(modify_st(pd$xlab), modify_st(pd$ylab)))
    }
    
    # Add estimate & SE
    coef.i$value <- pd$fit
    if (is.numeric(pd$se))
      coef.i$se  <- pd$se/pd$se.mult
    
    # Apply limits if requested
    tf.env <- environment(object$pcox$t.funcs[[
      which(sapply(object$pcox$smoothmap, function(x) i %in% x))
      ]])
    if (limit) {
      # Check if limits function exists
      limits.i <- tf.env$limits
      if (is.function(limits.i)) {
        # Apply limits function
        if (smooth.i$dim==2) {
          coef.i <- coef.i[limits.i(coef.i[[modify_st(pd$xlab)]],
                                    coef.i[[modify_st(pd$ylab)]]), ]
        } else if (smooth.i$dim==1) {
          coef.i <- coef.i[limits.i(coef.i[[modify_st(pd$xlab)]]), ]
        } else {
          stop("limit option only supported for 1-D and 2-D smooths")
        }
      }
    }
    coef.i
  })
  
  # Unlist if coef is length 1
  if(length(coef) == 1) {
    coef[[1]]
  } else {
    names(coef) <- sapply(object$smooth, function(x) x$label)
    coef
  }
}
# @rdname coefficients.pfr
#coef.pfr <- coefficients.pfr

