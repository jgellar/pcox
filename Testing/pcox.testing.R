# pcox.testing.R

# s(x)

N <- 500
J <- 101
x <- rnorm(N)
eta <- sin(x)
data1 <- simTVSurv(matrix(eta, nrow=N, ncol=J), Xdat=data.frame(x=x))
survobj <- Surv(data1$time, data1$event)

fit <- pcox(Surv(data1$time, data1$event) ~ s(data1$x))
fit <- pcox(Surv(time, event) ~ s(x), data=data1)
fit <- pcox(survobj ~ s(x), data=data1)

xidx <- seq(min(x), max(x), length=200)
bhat <- PredictMat(fit$pcox$smooth[[1]], data = data.frame(x=xidx)) %*% coef(fit)
plot(xidx, sin(xidx), type="l", lwd=3)
lines(xidx, bhat, col="red")
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
fit <- pcox(Surv(time,event) ~ lf(X1) + X2, data=data2)
bhat <- PredictMat(fit$pcox$smooth[[1]],
                   data = data.frame(X1.tmat=tind, L.X1=1)) %*% fit$coef[1:40]

# sofa data, with domain-standardized SOFA scores
data(sofa)
N <- nrow(sofa)
J <- ncol(sofa$SOFA)
sofa$SOFA.ds <- I(t(apply(sofa$SOFA, 1, function(x) {
  x.ind <- which(!is.na(x))
  if (length(x.ind)==1) {
    rep(x[x.ind],J)
  } else {
    approx(x.ind, x[x.ind], xout=seq(1, max(x.ind), length=J))$y
  }
})))

fit.sofa <- pcox(Surv())


X <- t(apply(X.mat, 1, function(x) {
  x.ind <- which(!is.na(x))
  approx(x.ind, x[x.ind], xout=seq(1, max(x.ind), length=J))$y
}))






# lf.vd





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


