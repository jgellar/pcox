#' @export

create.tt.func <- function(divide.by.t=FALSE, limits=NULL,
                           method=c("aic", "caic", "epic"),
                           integration=c("trapezoidal", "simpson", "riemann")) {
  method <- match.arg(method)
  integration <- match.arg(integration)
  
  if (is.null(limits))
    # Default
    limits <- "s<=t"
  if (!is.function(limits)) {
    if (!(limits %in% c("s<t", "s<=t")))
      stop("supplied <limits> argument unknown")
    if (limits == "s<t") {
      limits <- function(s, t) {s < t}
    } else {
      if (limits == "s<=t")
        limits <- function(s, t) {(s < t) | (s == t)}
    }
  }
  
  myfunc <- function(x.var,t.var,...) {
    n <- nrow(x.var)
    J <- ncol(x.var)
    tmat <- matrix(t.var, nrow=n, ncol=J)
    smat <- matrix(1:J, nrow=n, ncol=J, byrow=TRUE)
    mask <- t(outer(smat[1,], tmat[,1], limits))
    #L <- pcox:::getL(smat, integration=integration, n.int = tmat[,1])
    L <- getL3(smat, integration=integration, mask=mask)
    x.var[is.na(x.var)] <- 0
    rownames(x.var) <- NULL
    LX <- L*x.var
    if (divide.by.t) LX <- LX/tmat
    sm <- smoothCon(s(tmat, smat, by=LX),
                    data=data.frame(tmat=I(tmat), smat=I(smat), LX=I(LX)),
                    knots=NULL, absorb.cons=TRUE)[[1]]
    sm.out <<- sm
    pterm(sm, method=method, eps=.001)
  }
  #debug(myfunc)
  myfunc
}




