#' Extract coefficients from a fitted pcox model
#' 
#' 
#' 
#' 

coefficients.pcox <- function(object, raw=FALSE, term=NULL, inds=NULL, se=TRUE,
                              seWithMean=TRUE, n=100, n2=100, ...) {
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
      plotdata <- plotf(smooth.i, P=NULL, n=n, n2=n2, data=sdat.i, ...)
      
      if(is.re) plotdata$x <- levels(sdat.i[[smooth.i$term]])
      first <- smooth.i$first.para
      last  <- smooth.i$last.para
      coef.i <- list()
      coef.i$value <- drop(plotdata$X %*% object$coefficients[first:last])
      
      # ctrl-c-v from plot.mgcv.smooth :
      if (!is.null(plotdata$exclude)) coef.i$value[plotdata$exclude] <- NA
      if (se && plotdata$se) { ## get standard errors for fit
        ## test whether mean variability to be added to variability (only for centred terms)
        if (seWithMean && attr(smooth.i, "nCons") > 0) {
          
          ################# NEED TO TRANSLATE TO PCOX ############
          
          if (length(object$cmX) < ncol(object$Vp)){
            object$cmX <- c(object$cmX,rep(0,ncol(object$Vp)-length(object$cmX)))
          } 
          X1 <- matrix(object$cmX, nrow(plotdata$X), ncol(object$Vp), byrow=TRUE)
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
        if (!is.null(plotdata$exclude)) {
          coef.i$se[ plotdata$exclude] <- NA
          if (!is.null(object$var2)) coef.i$se2[plotdata$exclude] <- NA
        }
      }
      if (smooth.i$dim == 1) {
        if(!is.re) {
          coef.i[[gsub("\\.tmat", "\\.t", plotdata$xlab)]] <- plotdata$x
          #coef.i[[plotdata$xlab]] <- plotdata$x
        } else {
          coef.i[[smooth.i$term]] <- plotdata$x
        }
        
      } else {
        grid <- expand.grid(x=plotdata$x, y=plotdata$y)
        #coef.i[[plotdata$ylab]] <- grid$y
        #coef.i[[plotdata$xlab]] <- grid$x
        coef.i[[gsub("\\.tmat", "\\.t", plotdata$ylab)]] <- grid$y
        coef.i[[gsub("\\.smat", "\\.s", plotdata$xlab)]] <- grid$x
      }
      as.data.frame(coef.i)      
    })
    
    if (length(coefs)==1) coefs[[1]] else coefs
  } 
}

#' @rdname coefficients.pcox
coef.pcox <- coefficients.pcox