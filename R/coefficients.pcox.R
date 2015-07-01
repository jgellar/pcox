#' Extract coefficients from a fitted pcox model
#' 
#' This method is used to extract coefficients from a fitted `pcox` model, in
#' particular smooth functions resulting from including smooth effects
#' of scalar covariates, time-varying effects, functional predictors, or
#' historical effects.
#' 
#' @param object a fitted \code{pcox} object as produced by \code{\link{pcox}()}.
#' @param raw if \code{TRUE}, returns the parameters used in the fitting of the
#'   model. This would be the actual spline coefficients for smooth terms.
#' @param term integer to indicate the desired smooth term, according to its order
#'   in \code{object$pcox$smooth}. A vector of term indices may be supplied.
#' @param n number of estimation points for each 1-dimensional coefficient
#'   function.
#' @param n2 Square root of the number of points used to grid estimates of
#'   2-dimensinoal functions for contouring
#' @param se if \code{TRUE}, returns pointwise standard error estimates. Two
#'   estimates of the pointwise standard error are returned, based on the two
#'   estimates of the covariance matrix produced by \code{coxph}; see
#'   \code{\link[survival]{coxph}} for details.
#' @param seWithMean if \code{TRUE} the component smooths are shown with
#'   confidence intervals that include the uncertainty about the overall mean;
#'   if \code{FALSE}, then the uncertainty relates purely to the centered
#'   smooth itself.
#' @param exclude if \code{TRUE}, excludes reporting of the estimate at coordinates that are
#'   "too far" from data used to fit the model, as determined by
#'   \code{mgcv::plot.mgcv.smooth}, by setting the estimate to \code{NA}.
#' @param limit if \code{TRUE}, checks the \code{limits} function that was
#'   used to create the term and applies it to filter out some of the
#'   coordinates. Most relevant for variable-domain and historical functional
#'   terms.
#' @param ... other parameters to pass to the plotting function (either 
#'   \code{mgcv:::plot.mgcv.smooth} or \code{mgcv:::plot.random.effect}.
#' 
#' @return a data frame containing:
#'   \itemize{
#'     \item The coordinates where the coefficient is estimated
#'     \item The estimated coefficient at those coordinates (\code{value})
#'     \item The standard error estimate(s) at those coordinates, if requested
#'   }
#' 
#' @details
#'   This function calls either \code{mgcv:::plot.mgcv.smooth} or
#'   \code{mgcv:::plot.random.effect} to get the data that would be used
#'   to create an \code{mgcv}-style plot, but instead of plotting the data,
#'   it returns it. Much of the code and ideas behind this method were taken
#'   from \code{refund::coefficients.pfr}.
#' @author Jonathan Gellar <jgellar1@@jhu.edu> and Fabian Scheipl
#' @seealso \code{\link[mgcv]{plot.gam}}, \code{refund::coefficients.pfr}
#' @export
#' @importFrom stats coefficients

coefficients.pcox <- function(object, raw=FALSE, term=NULL, n=NULL, n2=NULL,
                              se=TRUE, seWithMean=TRUE, exclude=FALSE,
                              limit=TRUE, ...) {
  
  if (is.null(object$pcox$smooth) & !raw) {
    #warning("No smooth objects, returning raw coefficients")
    raw <- TRUE
  }
  if (raw) {
    # Raw coefficients only
    if (se) {
      ret <- list(coefficients=object$coefficients,
                  se  = sqrt(diag(object$var)))
      if (!is.null(object$var2)) ret$se2 <- sqrt(diag(object$var2))
      return(ret)
    } else {
      return(object$coefficients)
    }
  } else {
    # Return smooth/functional estimtes
    #term.labels <- attr(object$pcox$terms, "term.labels")
    if (is.null(term)) term <- seq_along(object$pcox$smooth)
    
    coefs <- lapply(term, function(i) {
      smooth.i <- object$pcox$smooth[[i]]
      sdat.i <- object$pcox$smoothdata[[i]]
      is.re <- "random.effect" %in% class(object$smooth.i)
      plotf <- if(is.re) {
        mgcv:::plot.random.effect
      } else {
        mgcv:::plot.mgcv.smooth
      }
      
      # Set values of n and n2 (if not supplied)
      n.i  <- ifelse(is.null(n),  ndefault(sdat.i[[smooth.i$term[[1]]]]), n)
      n2.i <- ifelse(is.null(n2), n.i, n2)
      
      # Get plotdata via plotf()
      plotdata <- plotf(smooth.i, P=NULL, n=n.i, n2=n2.i, data=sdat.i, ...)
      
      if(is.re) plotdata$x <- levels(sdat.i[[smooth.i$term]])
      first <- smooth.i$first.para
      last  <- smooth.i$last.para
      coef.i <- list()
      coef.i$value <- drop(plotdata$X %*% object$coefficients[first:last])
      
      # ctrl-c-v from plot.mgcv.smooth :
      if (exclude & !is.null(plotdata$exclude))
        coef.i$value[plotdata$exclude] <- NA
      if (se && plotdata$se) { ## get standard errors for fit
        ## test whether mean variability to be added to variability (only for centred terms)
        if (seWithMean && attr(smooth.i, "nCons") > 0) {
          if (length(object$means) < ncol(object$var)){
            object$means <- c(object$means, 
                              rep(0,ncol(object$var)-length(object$means)))
          }
          X1 <- matrix(object$means, nrow(plotdata$X), ncol(object$var),
                       byrow=TRUE)
          meanL1 <- smooth.i$meanL1
          if (!is.null(meanL1)) {
            X1 <- X1 / meanL1
          } 
          X1[,first:last] <- plotdata$X
          coef.i$se  <- sqrt(pmax(0, rowSums((X1 %*% object$var) * X1)))
          if (!is.null(object$var2))
            coef.i$se2 <- sqrt(pmax(0, rowSums((X1 %*% object$var2)* X1)))
        } else {
          ## se in centred (or anyway unconstained) space only
          coef.i$se  <- sqrt(pmax(0, rowSums((plotdata$X %*% 
            object$var[ first:last, first:last, drop=FALSE]) * plotdata$X)))
          if (!is.null(object$var2))
            coef.i$se2 <- sqrt(pmax(0, rowSums((plotdata$X %*% 
              object$var2[first:last, first:last, drop=FALSE]) * plotdata$X)))
        }
        if (exclude & !is.null(plotdata$exclude)) {
          coef.i$se[ plotdata$exclude] <- NA
          if (!is.null(object$var2)) coef.i$se2[plotdata$exclude] <- NA
        }
      }
      if (smooth.i$dim == 1) {
        if(!is.re) {
          coef.i[[modify_st(plotdata$xlab)]] <- plotdata$x
        } else {
          coef.i[[smooth.i$term]] <- plotdata$x
        }
        
      } else {
        grid <- expand.grid(x=plotdata$x, y=plotdata$y)
        coef.i[[modify_st(plotdata$ylab)]] <- grid$y
        coef.i[[modify_st(plotdata$xlab)]] <- grid$x
      }
      
      # Post-Processing
      coef.i <- as.data.frame(coef.i)
      tf.env <- environment(object$pcox$t.funcs[[
        which(sapply(object$pcox$smoothmap, function(x) i %in% x))
        ]])
      if (limit) {
        # Check if limits function exists
        limits.i <- tf.env$limits
        if (is.function(limits.i)) {
          # Apply limits function
          if (smooth.i$dim==2) {
            coef.i <- coef.i[limits.i(coef.i[[modify_st(plotdata$xlab)]],
                                      coef.i[[modify_st(plotdata$ylab)]]), ]
          } else if (smooth.i$dim==1) {
            coef.i <- coef.i[limits.i(coef.i[[modify_st(plotdata$xlab)]]), ]
          } else {
            stop("limit option only supported for 1-D and 2-D smooths")
          }
        }
      }
      coef.i
    })
    
    if (length(coefs)==1) coefs[[1]] else coefs
  } 
}

#' @rdname coefficients.pcox
coef.pcox <- coefficients.pcox

modify_st <- function(x) {
  ifelse(grepl("\\.tmat", x), "t", ifelse(grepl("\\.smat", x), "s", x))
}

ndefault <- function(x) {
  x <- unique(as.vector(x))
  x <- x[order(x)]
  if (length(unique(diff(x)))==1)
    length(x)
  else
    101
}
