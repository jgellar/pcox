#' @export

create.tt.func <- function(tind=NULL, divide.by.t=FALSE, limits=NULL,
                           method=c("aic", "caic", "epic"), eps=1e-6,
                           integration=c("trapezoidal", "simpson", "riemann"),
                           dbug=FALSE, k=NULL) {
  method <- match.arg(method)
  integration <- match.arg(integration)
  
  if (is.null(limits))
    # Default
    limits <- "s<=t"
  
  lag <- NULL
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
  sm.in <- NULL
  
  myfunc <- function(x.var,t.var,...) {
    n <- nrow(x.var)
    J <- ncol(x.var)
    if (is.null(tind))
      tind <- 1:J
    
    tmat <- matrix(t.var, nrow=n, ncol=J)
    smat <- matrix(tind, nrow=n, ncol=J, byrow=TRUE)
    mask <- t(outer(smat[1,], tmat[,1], limits))
    #L <- pcox:::getL(smat, integration=integration, n.int = tmat[,1])
    L <- pcox:::getL3(smat, integration=integration, mask=mask)
    x.var[is.na(x.var)] <- 0
    rownames(x.var) <- NULL
    LX <- L*x.var
    if (divide.by.t)
      LX <- if (is.null(lag)) {
        LX/matrix(t.var, nrow=n, ncol=J)
      } else {
        LX <- LX/matrix(sapply(t.var, min, lag), nrow=n, ncol=J)
      }
    pdat <- data.frame(tmat=I(tmat), smat=I(smat), LX=I(LX))
    
    if (is.null(sm.in)) {
      sm <- if (is.null(k)) {
        smoothCon(s(tmat, smat, by=LX), data=pdat,
                  knots=NULL, absorb.cons=TRUE)[[1]]
      } else {
        smoothCon(s(tmat, smat, by=LX, k=k), data=pdat,
                  knots=NULL, absorb.cons=TRUE)[[1]]
      }
      # NEED TO GET RID OF <<-
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


create.tt.smt.func2 <- function(divide.by.t=FALSE, limits=NULL,
                               method=c("aic", "caic", "epic"),
                               integration=c("trapezoidal", "simpson", "riemann"),
                               eps=1e-6, dbug=FALSE, sm.in=NULL, k=NULL,
                               rescale=F) {
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
    L <- pcox:::getL3(smat, integration=integration, mask=mask)
    LX <- L*x.var
    LX[!mask] <- NA
    
    J2 <- lag+1
    smtmat <- matrix((-lag:0), nrow=n, ncol=J2, byrow=TRUE)
    LX2 <- t(apply(LX, 1, function(x) {
      x <- x[!is.na(x)]
      c(rep(NA, J2-length(x)), x)
    }))
    LX2[is.na(LX2)] <- 0
    rownames(LX2) <- NULL
    tmat2 <- tmat[,1:J2]
    
    if (divide.by.t)
      LX2 <- LX2/matrix(sapply(t.var, min, lag), nrow=n, ncol=J2)
    
    if (rescale) {
      smtmat <- (smtmat-min(smtmat))/(max(smtmat)-min(smtmat))
      tmat2 <- (tmat2-min(tmat2))/(max(tmat2)-min(tmat2))
    }
    pdat <- data.frame(tmat=I(tmat2), smtmat=I(smtmat), LX=I(LX2))
    
    if (is.null(sm.in)) {
      sm <- if (is.null(k)) {
        smoothCon(s(tmat, smtmat, by=LX), data=pdat,
                  knots=NULL, absorb.cons=TRUE)[[1]]
      } else {
        smoothCon(s(tmat, smtmat, by=LX, k=k), data=pdat,
                  knots=NULL, absorb.cons=TRUE)[[1]]
      }
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




create.tt.smt.func <- function(divide.by.t=FALSE, limits=NULL,
                               method=c("aic", "caic", "epic"),
                               integration=c("trapezoidal", "simpson", "riemann"),
                               eps=1e-6, dbug=FALSE, sm.in=NULL, k=NULL,
                               rescale=F) {
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
    L <- pcox:::getL3(smat, integration=integration, mask=mask)
    x.var[is.na(x.var)] <- 0
    rownames(x.var) <- NULL
    LX <- L*x.var
    if (divide.by.t)
      LX <- LX/matrix(sapply(t.var, min, lag), nrow=n, ncol=J)
    
    smtmat <- smat-tmat
    if (rescale) {
      smtmat <- (smtmat-min(smtmat[mask]))/(max(smtmat[mask])-min(smtmat[mask]))
      tmat <- (tmat-min(tmat))/(max(tmat)-min(tmat))
    }
    pdat <- data.frame(tmat=I(tmat), smtmat=I(smtmat), LX=I(LX))
    
    if (is.null(sm.in)) {
      sm <- if (is.null(k)) {
        smoothCon(s(tmat, smtmat, by=LX), data=pdat,
                  knots=NULL, absorb.cons=TRUE)[[1]]
      } else {
        smoothCon(s(tmat, smtmat, by=LX, k=k), data=pdat,
                  knots=NULL, absorb.cons=TRUE)[[1]]
      }
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






create.tt.smtLine <- function(divide.by.t=FALSE, limits=NULL,
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
    
    # Intercept + slope
    cbind(rowSums(LX), rowSums(smtmat*LX))
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
    smtmat <- smat - tmat
    #L <- pcox:::getL(smat, integration=integration, n.int = tmat[,1])
    L <- pcox:::getL3(smat, integration=integration, mask=mask)
    x.var[is.na(x.var)] <- 0
    rownames(x.var) <- NULL
    LX <- L*x.var
    if (divide.by.t)
      LX <- LX/matrix(sapply(t.var, min, lag), nrow=n, ncol=J)
    pdat <- data.frame(smtmat=I(smtmat), LX=I(LX))
    
    if (is.null(sm.in)) {
      sm <- smoothCon(s(smtmat, by=LX), data=pdat,
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

