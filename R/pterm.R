if(getRversion() >= "2.15.1")  utils::globalVariables("theta")
#' Create a \code{coxph.penalty} object from an \code{mgcv}-style smooth
#' object
#' 
#' This function takes a an \code{mgcv}-style smooth object, created using
#' \code{\link[mgcv]{smoothCon}} or \code{\link[mgcv]{smooth.construct}},
#' and converts it into a \code{coxph.penalty} object that can be used
#' in a \code{survival::coxph} formula. It serves as the primary bridge
#' between the \code{mgcv} and \code{survival} packages.
#' 
#' @param sm the \code{mgcv}-style smooth object
#' @param method method for optimizing the smoothing parameter. Note that only
#'   \code{aic}, \code{caic}, and \code{epic} have been tested fully
#' @param eps tolerance level for optimizing the smoothing parameter
#' 
#' @return a \code{coxph.penalty} object, which is a design matrix with a
#'   number of atrributes that define the penalty. See Therneau and Grambsch
#'   (1988) for details.
#' @seealso \code{\link[survival]{coxph}} for a brief description of
#'   \code{coxph.penalty} terms; \code{\link[mgcv]{smoothCon}} and
#'   \code{\link[mgcv]{smooth.construct}} for describing \code{mgcv}-style
#'   smooth objects; \code{\link[survival]{pspline}},
#'   \code{\link[survival]{ridge}}, or \code{\link[survival]{frailty}} for
#'   other examples of \code{coxph.penalty} terms
#'   
#' @references
#' Therneau, T. M. and Grambsch, P. M. (1998). Penalized Cox models and Frailty.
#' Technical report, Division of Biostatistics. Mayo Clinic; Rochester, MN,
#' pages 1-58. Available at
#' http://www.mayo.edu/research/documents/frailtypdf/doc-10027273.
#' 
#' @author Jonathan Gellar <jgellar1@@jhu.edu>
#' @keywords internal
#' 

pterm <- function(sm, method=c("aic", "caic", "epic", "df", "fixed"), eps=1e-6) {
  method <- match.arg(method)
  W <- sm$X
  D <- sm$S[[1]]
  
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
                 df    = list(pfun=pfun.lFunc, cfun=frailty.controldf,
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

#' Control function for aic, caic, or epic-based optimization of a
#'   \code{coxph.penalty} term
#' 
#' @references
#' Therneau, T. M. and Grambsch, P. M. (1998). Penalized Cox models and Frailty.
#' Technical report, Division of Biostatistics. Mayo Clinic; Rochester, MN,
#' pages 1-58. Available at
#' http://www.mayo.edu/research/documents/frailtypdf/doc-10027273.
#' 
#' @keywords internal
control.aic <- function (parms, iter, old, n, df, loglik) {
  if (iter == 0) {
    if (is.null(parms$init)) 
      theta <- 0.005
    else theta <- parms$init[1]
    return(list(theta = theta, done = FALSE))
  }
  if (length(parms$type)) 
    type <- parms$type
  else type <- "aic"
  if (n < df + 2) 
    dfc <- (df - n) + (df + 1) * df/2 - 1
  else dfc <- -1 + (df + 1)/(1 - ((df + 2)/n))
  if (iter == 1) {
    history <- c(theta = old$theta, loglik = loglik, df = df, 
                 aic = loglik - df, aicc = loglik - dfc, epic = loglik-2*df)
    if (length(parms$init) < 2) 
      theta <- 1
    else theta <- parms$init[2]
    temp <- list(theta = theta, done = FALSE, history = history)
    return(temp)
  }
  history <- rbind(old$history, c(old$theta, loglik, df, loglik - 
                                    df, loglik - dfc, loglik - 2*df))
  if (is.null(parms$trace)) 
    trace <- FALSE
  else trace <- parms$trace
  if (iter == 2) {
    theta <- mean(history[, 1])
    return(list(theta = theta, done = FALSE, history = history, tst=4))
  }
  
  if (type=="caic") {
    aic <- history[,5]
  } else if (type=="epic") {
    aic <- history[,6]
  } else {
    aic <- history[,4]
  }
  # if (correct) 
  # aic <- history[, 5]
  # else aic <- history[, 4]
  done <- (abs(1 - aic[iter]/aic[iter - 1]) < parms$eps)
  x <- history[, 1]
  if (x[iter] == max(aic) && x[iter] == max(x)) 
    newtheta <- 2 * max(x)
  else newtheta <- frailty.brent(x, aic, lower = parms$lower, 
                                            upper = parms$upper)
  if (length(parms$trace) && parms$trace) {
    print(history)
    cat("    new theta=", format(newtheta), "\n\n")
  }
  list(theta = newtheta, done = done, history = history, tst=4)
}