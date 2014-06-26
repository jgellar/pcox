
lf.cox <- function(X, xind=NULL, fpca=FALSE, theta=NULL,
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
