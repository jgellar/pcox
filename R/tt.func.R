tt.func <- function(x.var,t.var,...) {
  n <- nrow(x.var)
  J <- ncol(x.var)
  tmat <- matrix(t.var, nrow=n, ncol=J)
  smat <- matrix(1:J, nrow=n, ncol=J, byrow=TRUE)
  L <- getL(smat, integration="trapezoidal", n.int=tmat[,1])
  x.var[is.na(x.var)] <- 0
  rownames(x.var) <- NULL
  LX <- L*x.var
  sm <- smoothCon(s(tmat, smat, by=LX), data=data.frame(tmat=I(tmat),
                                                        smat=I(smat),
                                                        LX=I(LX)),
                  knots=NULL, absorb.cons=TRUE)[[1]]
  pterm(sm, method="aic", eps=.001)
}

# fit <- coxph(Surv(Y,delta) ~ tt(X), tt=tt.func)
