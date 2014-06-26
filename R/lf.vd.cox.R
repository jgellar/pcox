lf.vd.cox <- function(X, tind = seq(0, 1, l = ncol(X)), Tind=NULL,
                  T.trans=identity,
                  domain=c("untransformed", "lagged", "standardized"),
                  interaction=c("nonparametric", "none", "linear", "quadratic"),
                  integration = c("simpson", "trapezoidal", "riemann"),
                  basistype=c("s","te","t2"),
                  rescale.unit = TRUE,
                  splinepars = NULL) {
  # splinepars = list(bs = "ps", k = min(ceiling(n/4), 40), m = c(2, 2))){
  
  n = nrow(X)
  nt = ncol(X)
  J.i <- apply(X, 1, function(x) max(which(!is.na(x))))
  domain <- match.arg(domain)
  interaction <- match.arg(interaction)
  integration <- match.arg(integration)
  basistype <- match.arg(basistype)
  tindname <- paste(deparse(substitute(X)), ".tmat", sep = "")
  Tindname <- paste(deparse(substitute(X)), ".Tmat", sep = "")
  LXname <- paste("L.", deparse(substitute(X)), sep = "")
  splinefun <- as.symbol(basistype)
  frmls <- if (is.null(splinepars)) {
    NULL
  } else {
    frmls <- formals(getFromNamespace(deparse(splinefun), ns = "mgcv"))
    modifyList(frmls[names(frmls) %in% names(splinepars)], 
               splinepars)
  }
  
  # Check domain/interaction/basis compatability
  if (domain %in% c("untransformed","lagged")) {
    if (interaction!="nonparametric") {
      stop("Untransformed and lagged domains require nonparametric interactions.")
    } else if (basistype %in% c("te","ti","t2")) {
      stop("Tensor product smooths are not supported for non-rectangular domains.")
    } else if (!is.null(splinepars)) {
      if (!is.null(splinepars$bs)) {
        if (!(splinepars$bs %in% c("ts","tp"))) {
          stop("Basis not supported for non-rectangular domains.")
        }
      }
    }
  } else if (interaction!="nonparametric" & basistype %in% c("te","ti","t2")) {
    stop("Tensor product smooths are not allowed for parametric interactions.")
  }
  
  # Create index matrices
  if (is.null(dim(tind))) {
    tind <- t(tind)
    stopifnot(ncol(tind) == nt)
    if (nrow(tind) == 1) {
      tind <- matrix(as.vector(tind), nrow = n, ncol = nt, 
                     byrow = T)
    }
    stopifnot(nrow(tind) == n)
  }
  if (is.null(Tind)) 
    Tind <- sapply(1:nrow(X), function(i) tind[i,max(which(!is.na(X[i,])))])
  Tind <- T.trans(Tind)
  if (is.null(dim(Tind))) {
    Tind <- t(Tind)
    stopifnot(ncol(Tind) == n)
    if (nrow(Tind) == 1) {
      Tind <- matrix(as.vector(Tind), nrow = n, ncol = nt)
    }
    stopifnot(nrow(tind) == n)
  }
  if (rescale.unit) {
    tind <- (tind-min(tind))/(max(tind)-min(tind))
    Tind <- (Tind-min(Tind))/(max(Tind)-min(Tind))
  }
  
  # Process functional predictor
  if (domain=="standardized") {
    X <- t(apply(X, 1, function(x) {
      J.i <- sum(!is.na(x))
      if (J.i==1) {
        rep(x[1], J)
      } else {
        approx(x=seq(0,1,length=J.i), y=x[1:J.i],
               xout=seq(0,1,length=J))$y
      }
    }))
    L <- getL(tind, integration=integration)
    LX <- L*X
  } else {
    L <- getL(tind, integration=integration, n.int=J.i)
    LX <- L*X
    if (domain=="lagged") {
      LX <- t(apply(LX, 1, function(x) {
        c(rep(NA,sum(is.na(x))), x[!is.na(x)])
      }))
    }
    LX[is.na(LX)] <- 0
  }
  
  if (interaction=="nonparametric") {
    data <- list(tind, Tind, LX)
    names(data) <- c(tindname, Tindname, LXname)
    call <- as.call(c(list(splinefun),
                      as.symbol(substitute(tindname)),
                      as.symbol(substitute(Tindname)),
                      by=as.symbol(substitute(LXname)), frmls))
  } else {
    data <- list(tind, LX)
    names(data) <- c(tindname, LXname)
    call <- as.call(c(list(splinefun), as.symbol(substitute(tindname)),
                      by=as.symbol(substitute(LXname)), frmls))		
    if (interaction %in% c("linear","quadratic")) {
      LX.lin <- LX * Tind
      LXname.lin <- paste0(LXname, ".lin")
      data[[LXname.lin]] <- LX.lin
      call <- call("+", call, as.call(c(list(splinefun), as.symbol(substitute(tindname)),
                                        by=as.symbol(substitute(LXname.lin)), frmls)))
    }
    if (interaction == "quadratic") {
      LX.qud <- LX * Tind^2
      LXname.qud <- paste0(LXname, ".qud")
      data[[LXname.qud]] <- LX.qud
      call <- call("+", call, as.call(c(list(splinefun), as.symbol(substitute(tindname)),
                                        by=as.symbol(substitute(LXname.qud)), frmls)))
    }
  }
  
  res <- list(call = call, data = data, L = L,
              tindname = tindname, Tindname=Tindname, LXname = LXname)
  return(res)
}






lf.vd.cox <- function(X, xind=NULL, fpca=FALSE, theta=NULL,
                      integration=c("trapezoidal","simpson","riemann"),
                      method=c("aic","caic","epic","reml", "df", "fixed"),
                      eps=.001, ...) {
  integration <- match.arg(integration)
  method <- match.arg(method)
  
  # Optional pre-smooth (most useful for misssing data)
  N <- nrow(X)
  X <- if (fpca) fpca.face.new(X)$Yhat else X
  X <- X - matrix(colMeans(X), nrow=N, ncol=ncol(X), byrow=TRUE)
  
  # Create design matrix W....
  if (is.null(xind)) {
    xind <- seq(0,1,length=ncol(X))
  }
  # sobj <- s(xind, ...)
  scon <- smoothCon(s(xind, ...), data=data.frame(xind=xind), knots=NULL, absorb.cons=FALSE)[[1]]
  phi.x <- scon$X
  D     <- scon$S[[1]]*scon$S.scale
  
  L.x <- getL(matrix(xind,nrow=N,ncol=length(xind),byrow=TRUE), integration)
  W <- (L.x * X) %*% phi.x
  
  # Penalty functions
  pfun.lFunc <- function(coef, theta, nevent, D) {
    lambda <- ifelse(theta<=0, 0, theta/(1-theta))
    list(penalty = as.numeric(t(coef) %*% D %*% coef) * lambda/2,
         first = lambda * D %*% coef,
         second = lambda * D, flag = FALSE)}
  
  # Control function based on optimization method
  temp <- switch(method, 
                 fixed = list(pfun=pfun.lFunc, cfun=function(parms, ...) {list(theta=parms$theta, done=TRUE)},
                              diag=FALSE, pparm=D, cparm=list(theta=theta)),
                 df    = list(pfun=pfun.lFunc, cfun=survival:::frailty.controldf,
                              
                              diag=FALSE, pparm=D),
                 reml  = list(pfun=pfun.lFunc, cfun=control.reml,
                              
                              diag=FALSE, pparm=D),
                 aic   = list(pfun=pfun.lFunc, cfun=control.aic,
                              cparm = list(eps=eps, init=c(0.5, 0.95), lower=0, upper=1, type="aic"),
                              diag=FALSE, pparm=D, cargs = c("neff", "df", "plik")),
                 caic  = list(pfun=pfun.lFunc, cfun=control.aic,
                              cparm = list(eps=eps, init=c(0.5, 0.95), lower=0, upper=1, type="caic"),
                              diag=FALSE, pparm=D, cargs = c("neff", "df", "plik")),
                 epic  = list(pfun=pfun.lFunc, cfun=control.aic,
                              cparm = list(eps=eps, init=c(0.5, 0.95), lower=0, upper=1, type="epic"),
                              diag=FALSE, pparm=D, cargs = c("neff", "df", "plik"))
  )
  
  # Return
  class(W) <- "coxph.penalty"
  attributes(W) <- c(attributes(W), temp)
  W
}
