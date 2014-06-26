
pp.cox <- function(W, D, theta=NULL, eps=.001,
                   method=c("aic","caic","epic","reml", "df", "fixed"), ...) {
  method <- match.arg(method)
  # model matrix W, penalty matrix D
  
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
