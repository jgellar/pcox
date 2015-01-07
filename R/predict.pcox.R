predict.pcox <- function(object, newdata=NULL, times=NULL,
                         type=c("lp", "risk", "expected"),
                         se.fit=FALSE, surv) {
  
  if (is.null(times)) {
    warning("times is NULL, returning predictions for training event times")
    times <- times.train
  }
  times <- unique(times)
  times <- times[order(times)]
  
  # STEP 1:  GET LINEAR PREDICTORS \eta_i(t)
  if (is.null(newdata)) {
    # Return predictions for training data
    if (all(times %in% unique(object$y[object$y[, 2]==1, 1]))) {
      # Already have linear predictors: need to put it in matrix format
      if (is.null(object$x)) stop("model must be fit with x=TRUE option")
      
      lp <- object$linear.predictors
      mapstrings <- strsplit(rownames(object$x), "[.]")
      subj <- as.numeric(sapply(mapstrings, function(x) x[1]))
      cols <- as.numeric(sapply(mapstrings, function(x) x[2]))
      cols <- ifelse(is.na(cols), 0, cols) + 1
      
      mapdat <- data.frame(subject=subj, col=cols, value=lp)
      lpmat  <- dcast(mapdat, subject~col)
      nc <- ncol(lpmat)-1
      lpmat[,2:(nc+1)] <- t(apply(lpmat[,-1], 1, function(x) {
        x0 <- x[!is.na(x)]
        c(rev(x0), rep(NA, nc-length(x0)))
      }))
      colnames(lpmat)[-1] <- rev(unique(object$y[,1]))
      
      #rset <- getRiskSet(surv)
      #rset$lp <- object$linear.predictors
    } else {
      # Need to calculate linear predictors at new times
      rset <- getRiskSet(surv, evaltimes = times)
      stop("Not yet supported")
      
      
      
      
    }
    
    
    
  } else {
    # New data supplied: calculate linear predictors for this data
    stop("Not yet supported!")
  }
}