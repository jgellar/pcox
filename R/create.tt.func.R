#' @export

create.tt.func <- function(divide=FALSE, limits=NULL) {
  if (is.null(limits)) limits <- function(s,t) {
    s<=t
  }
  function(x.var,t.var,...) {
    n <- nrow(x.var)
    J <- ncol(x.var)
    tmat <- matrix(t.var, nrow=n, ncol=J)
    smat <- matrix(1:J, nrow=n, ncol=J, byrow=TRUE)
    L <- getL(smat, integration="trapezoidal", n.int=tmat[,1])
    x.var[is.na(x.var)] <- 0
    rownames(x.var) <- NULL
    LX <- L*x.var
    if (divide) LX <- LX/tmat
    sm <- smoothCon(s(tmat, smat, by=LX), data=data.frame(tmat=I(tmat),
                                                          smat=I(smat),
                                                          LX=I(LX)),
                    knots=NULL, absorb.cons=TRUE)[[1]]
    sm.out <<- sm
    pterm(sm, method="aic", eps=.001)
  }
}
