#' Create a smooth and coxph.penalty object for use in a pcox formula
#' 
#' 

pcoxTerm <- function(data, limits, linear, tv, basistype, sind, 
                     integration, standardize, domain, basisargs,
                     method, eps, smooth, t=NULL) {
  
  # pcoxTerm will be called iff a smooth term is involved
  # data only includes the real (named) data (x)
  # time appears as t argument (vector of length nrow(data))
  # limits will be a function or NULL
  # if concurrent term, it will already be processed
  #   (so x_i(t) is in data as vector)
  
  # Returns either a list with the cpobj and smooth (if no input smooth),
  # or the prediction matrix (if there is an input smooth)
  
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
    
    
    # L <- getL4(sind, t, limits, standardize)
    
    #mask <- if (is.function(limits)) {
    #  t(outer(smat[1,], tmat[,1], limits))
    #} else NULL
    
    if (!is.null(t)) {
      tmat <- matrix(t, nrow=n, ncol=J)
      mask <- t(outer(smat[1,], tmat[,1], limits))
    } else {
      # t will be NULL if it's not a tt term: assume full range
      mask <- matrix(TRUE, nrow=n, ncol=J)
    }
    
    L <- getL3(smat, integration, mask)
    data[[1]][is.na(data[[1]])] <- 0
    evaldat$smat <- smat
    newcall <- c(newcall, quote(smat))
    if (tv) {
      # Allow smooth to vary over t
      evaldat$tmat <- tmat
      newcall <- c(newcall, quote(tmat))
    }
    if (!linear) {
      # Nonlinear functional term: add L matrix
      evaldat$L <- L
      newcall <- c(newcall, varnames, by=quote(L))
    } else if (length(data)==1) {
      # Linear functional term: add LX, remove X from smooth
      LX <- data[[1]]*L
      evaldat$LX <- LX
      newcall <- c(newcall, by=quote(LX))
    } else {
      stop("Not sure how to handle smooths with multiple linear functions yet")
    }
  }
  
  # Return is based on whether or not a smooth object has been supplied.
  # If not, we create the new smooth and the coxph.penalty object, and
  # return them. If smooth is supplied, return prediction matrix
  if (is.null(smooth)) {
    # Create and return smooth and cpobj objects
    newcall <- c(newcall, basisargs)
    smooth <- smoothCon(eval(as.call(newcall)), data=evaldat, knots=NULL,
                        absorb.cons=TRUE)
    if (length(smooth)>1) {
      # Can we turn these into a single smooth object (wider basis matrix, 
      #   block diagonal penalty matrix)?
      stop("We don't yet support terms with multiple smooth objects... stay tuned")
    }
    cpobj <- pterm(smooth[[1]], method=method, eps=eps)
    list(cpobj=cpobj, smooth=smooth)
  } else {
    # Return prediction matrix
    PredictMat(smooth[[1]], data=evaldat)
  }
}