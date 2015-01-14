predict.pcox <- function(object, type=c("lp", "risk", "expected"),
                         newdata=NULL, indices=NULL, ptimes=NULL, stimes=NULL,
                         n=NULL, se.fit=FALSE) {
  
  # newdata checks
  if (!is.null(newdata)) {
    if (!is.list(newdata))
      stop("newdata must be a list or data.frame")
    if (!is.data.frame(newdata)) {
      if (is.null(n)) {
        stop("n must be supplied if newdata is a list")
      }
    } else {
      if (is.null(n))
        n <- nrow(newdata)
      else if (n!=nrow(newdata)) {
        warning("n does not match extent of newdata - replacing n")
        n <- nrow(newdata)
      }
    }
  }
  
  any.tt <- any(object$pcox$t.types=="tt", na.rm = TRUE)
  lp.times <- rev(unique(object$y[,1]))
  
  # Calculate linear predictors
  lp <- if (is.null(newdata) &
              (!any.tt | is.null(ptimes) | all(ptimes %in% lp.times))) {
    # Predictions have already been calculated from training data
    if (!any.tt) {
      # No time-varying componenet: just return linear.predictor vector
      if (!is.null(ptimes))
        warning("ptimes ignored because model is not time-varying")
      object$linear.predictors
    } else {
      # Must organize linear.predictor vector into matrix
      if (is.null(ptimes)){
        warning(paste0("ptimes must be specified for time-varying model. ",
                       "Returning predictions for all training event times."))
        ptimes <- lp.times
      }
      if (is.null(object$x)) stop("model must be fit with x=TRUE option")
      
      # Generate map based on rownames of object$x
      lp <- object$linear.predictors
      mapstrings <- strsplit(rownames(object$x), "[.]")
      subj <- as.numeric(sapply(mapstrings, function(x) x[1]))
      cols <- as.numeric(sapply(mapstrings, function(x) x[2]))
      cols <- ifelse(is.na(cols), 0, cols) + 1
      
      # Convert to wide format, and put rows in correct order
      mapdat <- data.frame(subject=subj, col=cols, value=lp)
      lpmat  <- dcast(mapdat, subject~col)
      subj <- lpmat[,1]
      nc <- ncol(lpmat)-1
      lpmat <- t(apply(lpmat[,-1], 1, function(x) {
        x0 <- x[!is.na(x)]
        c(rev(x0), rep(NA, nc-length(x0)))
      }))
      colnames(lpmat) <- lp.times
      rownames(lpmat) <- subj
      
      # Filter to ptimes only
      if (length(ptimes) < length(lp.times))
        lpmat <- lpmat[, lp.times %in% ptimes]
      lpmat[, rank(ptimes)]
    }
  } else {
    # newdata or new ptimes are supplied - need to make new predictions
    if (is.null(newdata)) {
      newdata <- original.data # need this from pcox object?
    }
    if (any.tt) {
      # This is a time-varying model: calculate at specific ptimes
      # Need to expand dataset before predicting
      if (is.null(ptimes)) {
        # ptimes are needed: use training event times
        warning(paste0("ptimes must be specified for time-varying model. ",
                       "Returning predictions for all training event times."))
        ptimes <- lp.times
      }
      
      # Get survival times: needed to determine risk set at each of ptimes
      stimes <- if (is.null(stimes)) {
        if (is.null(newdata)) {
          # Use stimes from training data
          object$y[,1]
        } else {
          # Calculate stimes based on newdata
          stop("For now, need to enter stimes if using newdata")
          
          max.t <- max(ptimes)
          coltimes <- sapply(1:length(object$pcox$varmap), function(i) {
            vnames <- object$pcox$varmap[[i]]
            
          #coltimes <- sapply(object$pcox$varmap, function(vnames) {
            sapply(vnames, function(vn.i) {
              var.i <- newdata[[vn.i]]
              if (is.matrix(var.i)) {
                idx.i <- if (!is.null(indices[[vn.i]]))
                  vn.i
                else something # get from tt function? ugh
              }
            })
            
            vars <- newdata[[vnames]]
            maxt <- sapply(vars, function(v) {
              if (is.matrix(v)) {
                idx.v <- 
                apply(v, 1, function(x) max(which(!is.na(x))))
              }
            })
            # do something where you go through vars, using indices, and calculate
            # max(which(!is.na(x))) for each matrix
          })
        }
      } else {
        # stimes is available - return (depending on type)
        if (class(stimes)=="Surv")
          stimes[,1]
        else if (is.numeric(stimes))
          stimes
        else
          stop("Unrecognized type for stimes")
      }
      
      if (length(stimes)!=n) {
        stop("length of stimes does not match n")
      }
      
      # Expand dataset at ptimes
      expand.map <- do.call("rbind", sapply(ptimes, function(t.i) {
        cbind(t.i, which(stimes>=t.i))
      }))
      ptimes <- expand.map[,1]
      #newdata <- newdata[expand.map[,2]]
    }
    
    t.funcs <- object$pcox$t.funcs
    t.types <- object$pcox$t.types
    
    # Get prediction matrix (X for \hat\eta = X %*% \hat\beta)
    pmat <- do.call("cbind", lapply(1:length(t.funcs), function(i) {
      vnames <- object$pcox$varmap[[i]]
      
      if (is.na(t.types[i])) {
        # No transformation: check if data is available, and if so, return
        if (!all(vnames %in% names(newdata))) {
          v.miss <- vnames[!(vnames %in% names(newdata))]
          stop(paste0("Variable", ifelse(length(v.miss)==1, " ", "s "),
                      paste(v.miss, collapse = ", "),
                      ifelse(length(v.miss)==1, " is ", " are "),
                      "required but not supplied"))
        }
        p.i <- as.matrix(newdata[vnames])
        if (any.tt) p.i[expand.map[,2], ] else p.i
      } else {
        # Check if newdata is supplied using original varialbes or in mgcv form
        mnames <- strsplit(object$pcox$labelmap[[i]], "[/(]|,|[/)]:")[[1]][-1]
        if (all(vnames %in% names(newdata))) {
          # newdata is supplied using original variable names/format
          # use t.func
          x.i <- newdata[vnames]
          switch(t.types[i], xt = {
            # Execute xt-function, then possibly expand result
            p.i <- t.funcs[[i]](x.i)
            if (any.tt) p.i[expand.map[,2], ] else p.i
          }, tt = {
            # Execute tt-function on expanded data
            t.funcs[[i]](x.i[expand.map[,2], ], ptimes)
          }, stop(paste0("Unrecognized t.type for term ", i)))
        } else if (all(mnames %in% names(newdata))) {
          # newdata is supplied using mgcv-style format
          # Call PredictMat() directly
          PredictMat(environment(t.funcs[[i]])$smooth[[1]], newdata[mnames], n)
        } else {
          stop(paste0("Required data is not suppled for term ", i))
        }
      }
    }))
    
    # Calculate X %*% \hat\beta
    lp <- pmat %*% object$coefficients - 
      as.numeric(object$means %*% object$coefficients)
    
    if (any.tt) {
      # Put in matrix format based on expand.map
      acast(data.frame(ptimes=ptimes, subject=expand.data[,2], lp=lp),
            subject ~ ptimes, value.var=lp)
    } else as.vector(lp)
  }
  
  lp
}