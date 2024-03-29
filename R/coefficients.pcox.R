#' Extract coefficient functions from a fitted pcox-object
#' 
#' This function is used to extract a coefficient from a fitted `pcox` model, in
#' particular smooth functions resulting from including smooth effects
#' of scalar covariates, time-varying effects, functional predictors, or
#' historical effects. Only one coefficient can be extracted at a time.
#' 
#' @param object return object from \code{\link{pcox}}
#' @param select integer indicating the index of the desired smooth term
#'   in \code{object$pcox$smooth}. Enter 0 to request the raw coefficients
#'   (i.e., \code{object$coefficients}) and standard errors (if \code{se==TRUE}).
#' @param coords named list indicating the desired coordinates where the
#'   coefficient function is to be evaluated. Names must match the names of
#'   \code{object$pcox$smooth[[selet]]$term}. If \code{NULL}, uses \code{n}
#'   to generate equally-spaced coordinates.
#' @param n integer vector indicating the number of equally spaced coordinates
#'   for each term. If length 1, the same number is used for each term.
#'   Otherwise, the length must match \code{object$pcox$smooth[[select]]$dim}.
#' @param se if \code{TRUE}, returns pointwise standard error estimates. Defaults
#'   to \code{FALSE} if raw coefficients are being returned; otherwise \code{TRUE}.
#' @param seWithMean if \code{TRUE} the standard errors include uncertainty about
#'   the overall mean; if \code{FALSE}, they relate purely to the centered
#'   smooth itself. Marra and Wood (2012) suggests that \code{TRUE} results in
#'   better coverage performance for GAMs, though we have not tested this for
#'   survival models.
#' @param limit if \code{TRUE}, checks if a \code{limits} function was
#'   used to generate the term, and if so, applies the function to the
#'   output to only produce estimates over regions that were within "limits".
#' @param ... these arguments are ignored
#'   
#' @details
#' Should add a description here about entering coordinates.
#'   
#' @return a data frame containing the evaluation points, 
#'    coefficient function values and optionally the SE's for the term indicated
#'    by \code{select}.
#' 
#' @author Fabian Scheipl and Jonathan Gellar
#' 
#' @references
#' Marra, G and S.N. Wood (2012) Coverage Properties of Confidence Intervals for
#' Generalized Additive Model Components. Scandinavian Journal of Statistics.
#' 
#' @importFrom stats coefficients
#' @export

coefficients.pcox <- function(object, select=1, coords=NULL, n=NULL,
                   se=ifelse(length(object$pcox$smooth) & select, TRUE, FALSE),
                   seWithMean=FALSE, limit=TRUE, ...) {
  
  if (!length(object$pcox$smooth) | select==0) {
    if (se) {
      ret <- list(coefficients=object$coefficients,
                  se  = sqrt(diag(object$var)))
      if (!is.null(object$var2)) ret$se2 <- sqrt(diag(object$var2))
      return(ret)
    } else {
      return(object$coefficients)
    }
  }
  
  object$coefficients[is.na(object$coefficients)] <- 0 #Zero out singular coefs
  smooth.i <- object$pcox$smooth[[select]]
  sdat.i   <- object$pcox$smoothdata[[select]]
  
  coef.i <- if ("random.effect" %in% class(smooth.i)) {
    stop("Random effects not yet implemented for coef.pcox")
    
  } else {
    # Check or create coords
    if (is.null(coords)) {
      if (is.null(n))
        n <- sapply(sdat.i[smooth.i$term], ndefault)
      else if (length(n)==1)
        n <- rep(n, smooth.i$dim)
      else if (length(n)!=smooth.i$dim)
        stop("length of n must match the number of terms of the smooth")
      coords <- mapply(function(x,y) {
        seq(min(x), max(x), length=y)
      }, sdat.i[smooth.i$term], n, SIMPLIFY=FALSE)
    } else if (!all(names(coords) %in% smooth.i$term)) {
      stop("coords must be a list with names equal to smooth[[select]]$term")
    }
    
    # Create grid and get prediction matrix
    coef.i <- do.call(expand.grid, coords)[smooth.i$term]
    if (smooth.i$by!="NA")
      coef.i[smooth.i$by] <- 1
    pmat <- mgcv::PredictMat(smooth.i, coef.i)
    
    # Get coefficient function
    coef.i[smooth.i$by] <- NULL
    first <- smooth.i$first.para
    last  <- smooth.i$last.para
    coef.i$value <- pmat %*% object$coefficients[first:last]
    
    #### GET SE: COPIED FROM plot.mgcv.smooth
    if (se) { ## get standard errors for fit
      ## test whether mean variability to be added to variability (only for centred terms)
      if (seWithMean && attr(smooth.i, "nCons")>0) {
        if (length(object$means) < ncol(object$var)){
          object$means <- c(object$means, 
                            rep(0,ncol(object$var)-length(object$means)))
        }
        
        X1 <- matrix(object$means, nrow(coef.i), ncol(object$var), byrow=TRUE)
        meanL1 <- smooth.i$meanL1
        if (!is.null(meanL1)) {
          X1 <- X1 / meanL1
        } 
        X1[,first:last] <- pmat
        coef.i$se  <- sqrt(pmax(0, rowSums((X1 %*% object$var) * X1)))
        if (!is.null(object$var2))
          coef.i$se2 <- sqrt(pmax(0, rowSums((X1 %*% object$var2)* X1)))
        
      } else {
        ## se in centred (or anyway unconstained) space only
        coef.i$se  <-
          sqrt(pmax(0, rowSums((pmat %*% object$var[ first:last, first:last,
                                                     drop=FALSE]) * pmat)))
        if (!is.null(object$var2))
          coef.i$se2 <-
          sqrt(pmax(0, rowSums((pmat %*% object$var2[first:last, first:last,
                                                     drop=FALSE]) * pmat)))
      }
      
    }
    
    # Apply limits if requested
    if (limit) {
      # Check if limits function exists
      tf.env <- environment(object$pcox$t.funcs[[
        which(sapply(object$pcox$smoothmap, function(x) select %in% x))
        ]])
      limits.i <- tf.env$limits
      if (is.function(limits.i)) {
        args <- coef.i[1:length(formals(limits.i))]
        names(args) <- names(formals(limits.i))
        coef.i <- coef.i[do.call(limits.i, args), ]
      }
    }
    coef.i
  }
  
  names(coef.i)[1:smooth.i$dim] <- sapply(smooth.i$term, modify_st)
  coef.i
}


#' @rdname coefficients.pcox
#' @export
coef.pcox <- coefficients.pcox

modify_st <- function(x) {
  ifelse(grepl("\\.tmat", x), "t", ifelse(grepl("\\.smat", x), "s", x))
}

ndefault <- function(x) {
  x <- unique(as.vector(x))
  x <- x[order(x)]
  
  diffs <- unique(round(diff(x), 10))
  n <- if (length(diffs)<5) {
    length(seq(min(x), max(x), by=min(diffs)))
  } else 101
  while(n>101)
    n <- round(n/2)
  n
}
