library(devtools)
library(mgcv)
library(survival)

install.packages("/Users/jonathangellar/Documents/Projects/refundRepos/refund/pkg/refundDevel", repos=NULL, type="source")
library(refundDevel)


dev_mode()
load_all()


# pcox.testing.R

# s(x)

N <- 500
J <- 101
x <- rnorm(N)
eta <- sin(x)
data1 <- simTVSurv(matrix(eta, nrow=N, ncol=J), Xdat=data.frame(x=x))
survobj <- Surv(data1$time, data1$event)

fit <- pcox(Surv(time, event) ~ s(x), data=data1, eps = .00000001)
fit1<- pcox(Surv(time, event) ~ s(x), data=data1, eps = .00000001)

fit <- pcox(survobj ~ s(x), data=data1)

xidx <- seq(min(x), max(x), length=200)
bhat <- PredictMat(fit$pcox$smooth[[1]], data = data.frame(x=xidx)) %*% coef(fit)
bhat1<- PredictMat(fit1$pcox$smooth[[1]], data = data.frame(x=xidx)) %*% coef(fit1)


plot(xidx, sin(xidx), type="l", lwd=3)
lines(xidx, bhat, col="red")
lines(xidx, bhat1,col="blue")
rug(x)


# lf(X)
X <- genX(N)
X <- X - matrix(colMeans(X), nrow=N, ncol=J, byrow=TRUE)
X2 <- rnorm(N)
Xdat <- data.frame(X1=I(X), X2=X2)
class(Xdat$X1) <- NULL
eta <- apply(Xdat, 1, function(x) {
  mean(beta * x[1:J]) + x[J+1]
})
data2 <- cbind(simTVSurv(matrix(eta, nrow=N, ncol=J)), Xdat)

fit <- pcox(Surv(time,event) ~ lf(X1) + X2, data=data2, eps=1e-6)
fit1<- pcox(Surv(time,event) ~ lf(X1) + X2, data=data2, eps=1e-6)

bhat <- PredictMat(fit$pcox$smooth[[1]],
                   data = data.frame(X1.tmat=tind, L.X1=1)) %*% fit$coef[1:40]
bhat1<- PredictMat(fit1$pcox$smooth[[1]],
                   data = data.frame(X1.tmat=tind, L.X1=1)) %*% fit1$coef[1:40]

plot(tind, beta, type="l", lwd=3)
lines(tind, bhat, col="red")
lines(tind, bhat1, col="blue")

# sofa follow-up data, with domain-standardized SOFA scores
data(sofa_fu)
N <- nrow(sofa_fu)
J <- ncol(sofa_fu$SOFA)
class(sofa_fu[[3]]) <- NULL
class(sofa_fu[[4]]) <- NULL
class(sofa_fu[[5]]) <- NULL
sofa_fu$log.los <- log(sofa_fu$los)

xind <- seq(log(2), log(157), length=101)
tind <- seq(0,1,length=101)
fhat.dat <- data.frame(log.los=xind)
bhat.dat <- data.frame(SOFA_ds.tmat=tind, L.SOFA_ds=1)

fits.sofa <- list()
fhat = bhat <- matrix(nrow=3, ncol=101)
Zhat <- matrix(nrow=3, ncol=3)
methods <- c("aic", "caic", "epic")
for (m in 1:3) {
  fit <- pcox(Surv(time,event) ~ age + male + Charlson + s(log.los)
                   + lf(SOFA_ds, splinepars=list(bs="ps", m=c(2,2))),
                   data=sofa_fu, method=methods[m], eps=.000001)
  fits.sofa[[m]] <- fit
  fhat[m,] <- PredictMat(fit$pcox$smooth[[1]], data=fhat.dat) %*% fit$coef[4:12]
  bhat[m,] <- PredictMat(fit$pcox$smooth[[2]], data=bhat.dat) %*% fit$coef[13:52]
  Zhat[m,] <- fit$coef[1:3]
}


B <- 10
samp.boot <- matrix(nrow=B, ncol=N)
fhat.boot <- array(dim=c(B,3,length(xind)))
bhat.boot <- array(dim=c(B,3,length(tind)))
Zhat.boot <- array(dim=c(B,3,3))

for (b in 1:B) {
  print(b)
  samp <- sample(1:N, replace=TRUE)
  samp.boot[b,] <- samp
  dat.b <- sofa_fu[samp,]
  for (m in 1:3) {
    fit.b <- pcox(Surv(time,event) ~ age + male + Charlson + s(log.los)
                  + lf(SOFA_ds), data=dat.b, method=methods[m],
                  eps=.0000001)
    fhat.boot[b,m,] <- PredictMat(fit.b$pcox$smooth[[1]], data=fhat.dat
    ) %*% fit.b$coef[4:12]
    bhat.boot[b,m,] <- PredictMat(fit.b$pcox$smooth[[2]], data=bhat.dat
    ) %*% fit.b$coef[13:52]
    Zhat.boot[b,m,] <- fit.b$coef[1:3]
  }
  
}


estplot <- function(i, m) {
  if (i==1) {
    est <- fhat[m,]
    boot<- fhat.boot[,m,]
    inds<- xind
    ylm <- NULL
  } else {
    est <- bhat[m,]
    boot<- bhat.boot[,m,]
    inds<- tind
    ylm <- c(-10,10)
  }
  plot(inds, est, type="l", lwd=3, ylim=ylm)
  for (b in 1:nrow(boot)) {
    lines(inds, boot[b,], col=b, lwd=.5)
  }
  lines(inds, est, lwd=3)
}










# lf.vd

fit <- pcox(Surv(time,event) ~ age + male + Charlson + s(log.los)
            + lf.vd(SOFA, splinepars=list(k=50)),
            data=sofa_fu, method=methods[m], eps=.000001)
LX.pre  <- matrix(nrow=J,ncol=J)
LX.pre[lower.tri(LX.pre,diag=TRUE)] <- 1
pre.chk <- data.frame(SOFA.tmat=I(matrix(seq(0,1,length=J), nrow=J, ncol=J, byrow=TRUE)),
                      SOFA.Tmat=I(matrix(seq(0,1,length=J), nrow=J, ncol=J)),
                      L.SOFA=I(LX.pre))
fit$smooth <- fit$pcox$smooth
est1 <- getEst(fit, pre.chk, sm.id = 2)

est <- getEst(fit, data=)


fits.sofa.vd <- list()
fhat.vd = bhat.vd <- matrix(nrow=3, ncol=101)
Zhat.vd <- matrix(nrow=3, ncol=3)
methods <- c("aic", "caic", "epic")
for (m in 1:3) {
  fit <- pcox(Surv(time,event) ~ age + male + Charlson + s(log.los)
              + lf.vd(SOFA_ds, splinepars=list(bs="ps", m=c(2,2))),
              data=sofa_fu, method=methods[m], eps=.000001)
  fits.sofa.vd[[m]] <- fit
  fhat.vd[m,] <- PredictMat(fit$pcox$smooth[[1]], data=fhat.dat) %*% fit$coef[4:12]
  bhat.vd[m,] <- PredictMat(fit$pcox$smooth[[2]], data=bhat.dat) %*% fit$coef[13:52]
  Zhat.vd[m,] <- fit$coef[1:3]
}



T.r <- -log(runif(N)) * exp(-eta)
Y <- (w.scale^.75 * T.r)^(1/.75)
C <- rep(365*2, N)
delta <- Y<C
Y[!delta] <- C[!delta]
data2 <- data.frame(time=Y, event=delta)
data2 <- cbind(data2, Xdat)





Xdat <- data.frame(X1=I(genX(N)), X2=rnorm(N))
tind <- seq(0,1,length=101)
beta <- sin(2*pi*tind)






eta  <- t(apply(Xdat, 1, function(x) {
  0.5*x[-length(x)] + x[length(x)]
}))
data2 <- simTVSurv(eta, Xdat)


