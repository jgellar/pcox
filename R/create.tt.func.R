#' @export

create.tt.func <- function(divide.by.t=FALSE, limits=NULL,
                           method=c("aic", "caic", "epic"),
                           integration=c("trapezoidal", "simpson", "riemann"),
                           eps=1e-6, dbug=FALSE, sm.in=NULL) {
  method <- match.arg(method)
  integration <- match.arg(integration)
  
  if (is.null(limits))
    # Default
    limits <- "s<=t"
  
  if (!is.function(limits)) {
    if (is.numeric(limits)) {
      # limits specifies a lag for inclusion
      lag <- limits
      limits <- function(s, t) {s<=t & s>=(t-lag)}
    } else if (limits == "s<t") {
      limits <- function(s, t) {s < t}
    } else if (limits == "s<=t") {
      limits <- function(s, t) {(s <= t)}
    } else {
      stop("supplied <limits> argument unknown")
    }
  }
  
  myfunc <- function(x.var,t.var,...) {
    n <- nrow(x.var)
    J <- ncol(x.var)
    tmat <- matrix(t.var, nrow=n, ncol=J)
    smat <- matrix(1:J, nrow=n, ncol=J, byrow=TRUE)
    mask <- t(outer(smat[1,], tmat[,1], limits))
    #L <- pcox:::getL(smat, integration=integration, n.int = tmat[,1])
    L <- pcox:::getL3(smat, integration=integration, mask=mask)
    x.var[is.na(x.var)] <- 0
    rownames(x.var) <- NULL
    LX <- L*x.var
    if (divide.by.t)
      LX <- LX/matrix(sapply(t.var, min, lag), nrow=n, ncol=J)
    pdat <- data.frame(tmat=I(tmat), smat=I(smat), LX=I(LX))
    
    if (is.null(sm.in)) {
      sm <- smoothCon(s(tmat, smat, by=LX), data=pdat,
                      knots=NULL, absorb.cons=TRUE)[[1]]
      sm.out <<- sm
      pterm(sm, method=method, eps=eps)
    } else {
      PredictMat(sm.in, data = pdat)
    }
  }
  if (dbug) {
    debug(myfunc)
  } else if (isdebugged(myfunc))
    undebug(myfunc)
  myfunc
}


create.tt.smtlinear <- function(divide.by.t=FALSE, limits=NULL,
                                method=c("aic", "caic", "epic"),
                                integration=c("trapezoidal", "simpson", "riemann"),
                                eps=1e-6, dbug=FALSE) {
  method <- match.arg(method)
  integration <- match.arg(integration)
  lag <- limits
  limits <- function(s,t) {s<=t & s>=(t-lag)}
  
  myfunc <- function(x.var,t.var,...) {
    n <- nrow(x.var)
    J <- ncol(x.var)
    tmat <- matrix(t.var, nrow=n, ncol=J)
    smat <- matrix(1:J, nrow=n, ncol=J, byrow=TRUE)
    mask <- t(outer(smat[1,], tmat[,1], limits))
    #L <- pcox:::getL(smat, integration=integration, n.int = tmat[,1])
    L <- pcox:::getL3(smat, integration=integration, mask=mask)
    x.var[is.na(x.var)] <- 0
    rownames(x.var) <- NULL
    
    smtmat <- smat-tmat
    LX <- L*x.var
    if (divide.by.t)
      LX <- LX/matrix(sapply(t.var, min, lag), nrow=n, ncol=J)
    rowSums(smtmat*LX)
  }
  if (dbug) {
    debug(myfunc)
  } else if (isdebugged(myfunc))
    undebug(myfunc)
  myfunc
}



create.tt.func.ntv <- function(divide.by.t=FALSE, limits=NULL,
                           method=c("aic", "caic", "epic"),
                           integration=c("trapezoidal", "simpson", "riemann"),
                           eps=1e-6, dbug=FALSE, sm.in=NULL) {
  method <- match.arg(method)
  integration <- match.arg(integration)
  
  if (is.null(limits))
    # Default
    limits <- "s<=t"
  
  if (!is.function(limits)) {
    if (is.numeric(limits)) {
      # limits specifies a lag for inclusion
      lag <- limits
      limits <- function(s, t) {s<=t & s>=(t-lag)}
    } else if (limits == "s<t") {
      limits <- function(s, t) {s < t}
    } else if (limits == "s<=t") {
      limits <- function(s, t) {(s <= t)}
    } else {
      stop("supplied <limits> argument unknown")
    }
  }
  
  myfunc <- function(x.var,t.var,...) {
    n <- nrow(x.var)
    J <- ncol(x.var)
    tmat <- matrix(t.var, nrow=n, ncol=J)
    smat <- matrix(1:J, nrow=n, ncol=J, byrow=TRUE)
    mask <- t(outer(smat[1,], tmat[,1], limits))
    #L <- pcox:::getL(smat, integration=integration, n.int = tmat[,1])
    L <- pcox:::getL3(smat, integration=integration, mask=mask)
    x.var[is.na(x.var)] <- 0
    rownames(x.var) <- NULL
    LX <- L*x.var
    if (divide.by.t)
      LX <- LX/matrix(sapply(t.var, min, lag), nrow=n, ncol=J)
    pdat <- data.frame(smat=I(smat), LX=I(LX))
    
    if (is.null(sm.in)) {
      sm <- smoothCon(s(smat, by=LX), data=pdat,
                      knots=NULL, absorb.cons=TRUE)[[1]]
      sm.out <<- sm
      pterm(sm, method=method, eps=eps)
    } else {
      PredictMat(sm.in, data = pdat)
    }
  }
  if (dbug) {
    debug(myfunc)
  } else if (isdebugged(myfunc))
    undebug(myfunc)
  myfunc
}

