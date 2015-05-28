#' Extract coefficients from a fitted pcox model
#' 
#' 
#' 
#' 

coefficients.pcox <- function(object, raw=FALSE, term=NULL, n=100, n2=100,
                              inds=NULL, se=TRUE,
                              limit=TRUE, untransform=TRUE,
                              seWithMean=TRUE, ...) {
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
          coef.i[[modify_st(plotdata$xlab)]] <- plotdata$x
          #coef.i[[gsub("\\.tmat", "\\.t", plotdata$xlab)]] <- plotdata$x
          #coef.i[[plotdata$xlab]] <- plotdata$x
        } else {
          coef.i[[smooth.i$term]] <- plotdata$x
        }
        
      } else {
        grid <- expand.grid(x=plotdata$x, y=plotdata$y)
        coef.i[[modify_st(plotdata$ylab)]] <- grid$y
        coef.i[[modify_st(plotdata$xlab)]] <- grid$x
        #coef.i[[plotdata$ylab]] <- grid$y
        #coef.i[[plotdata$xlab]] <- grid$x
        #coef.i[[gsub("\\.tmat", "\\.t", plotdata$ylab)]] <- grid$y
        #coef.i[[gsub("\\.smat", "\\.s", plotdata$xlab)]] <- grid$x
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
          if (smooth.i$dim==2 & "s" %in% names(coef.i) & "t" %in% names(coef.i))
            coef.i <- coef.i[limits.i(coef.i$s, coef.i$t), ]
        } else {
          stop("limits option not yet supported for this term")
        }
        
      }
      if (untransform) {
        # Check if transformation functions were applied
        st.i <- tf.env$s.transform
        tt.i <- tf.env$t.transform
        if (!is.null(st.i)) {
          stop("untransform option not yet supported")
        }
        if (!is.null(tt.i)) {
          stop("untransform option not yet supported")
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
