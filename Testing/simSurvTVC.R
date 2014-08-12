library(timereg)
library(mgcv)

# Generate random survival data based on a PH model with
# TVCs or TVEs

# SCENARIO 1: time-varying effects of scalar covariates
N <- 500
Xdat <- data.frame(X1=rbinom(N,1,.5), X2=rnorm(N))
tms  <- seq(0, 2*pi, length=100)
beta <- 3*sin(tms)
eta  <- t(apply(Xdat, 1, function(x) {
  -2 + 4*x[1] + beta*x[2]
}))
data1 <- simTVSurv(eta, Xdat)

# Model 1:  timecox
fit1.tc <- timecox(Surv(time,event) ~ const(X1) + X2, data1)
tmp.t <- fit1.tc$cum[,1]
tmp.x <- fit1.tc$cum[,3]
est1.tc <- diff(predict(gam(tmp.x~s(tmp.t))))/diff(tmp.t)
plot(est1.tc ~ tmp.t[-length(tmp.t)], type="l", ylim=c(-4,4))
lines(beta ~ c(1:100), col="red")

# Model 2:  coxph+tt()
fit1.tt <- coxph(Surv(time,event) ~ X1 + tt(X2), data=data1,
                 tt=function(x,t,...) x*s.cox(t))
sm <- smoothCon(s(tvec), data=data.frame(tvec=fit1.tt$y[,1]),
                knots=NULL, absorb.cons=TRUE)[[1]]
pmat <- PredictMat(sm, data=data.frame(tvec=1:100))
est1.tt <- as.vector(pmat %*% coef(fit1.tt)[-1])
lines(as.vector(est1.tt) ~ c(1:100), col="blue")


# SCENARIO 2: Time-Varying Covariates - concurrent effect
Xdat <- data.frame(X1=I(genX(N)), X2=rnorm(N))
eta  <- t(apply(Xdat, 1, function(x) {
  0.5*x[-length(x)] + x[length(x)]
}))
data2 <- simTVSurv(eta, Xdat)
fit2 <- coxph(Surv(time,event) ~ tt(X1) + X2, data=data2,
              na.action=na.pass,
              tt=function(x,t,...) {sapply(1:length(t), function(i)
                x[i,t[i]])
                })


# SCENARIO 3: Historical Functional Terms for TVC's
N <- 500
J <- 101
Xdat <- data.frame(X1=I(genX(N, s=seq(0,1,length=J))), X2=rnorm(N))
beta1 <- makeBetaMat(J, genBeta1)
eta <- sapply(1:J, function(j) {
  1.3*Xdat[[2]] + Xdat[[1]][,1:j,drop=F] %*% beta1[j,1:j] / j
})
data3 <- simTVSurv(eta, Xdat=Xdat)
fit3 <- coxph(Surv(time,event) ~ tt(X1) + X2, data=data3, na.action=na.pass,
              tt=tt.funcD)
est3 <- getHCEst(sm.out, 1:101, coefs = coef(fit3)[-31])

par(mfrow=c(1,2))
image(t(beta1), zlim=c(-6,6), col=jet.colors(64), xaxt="n", yaxt="n", main="True Beta")
axis(1, seq(0,1,length=5), seq(0,100,length=5))
axis(2, seq(0,1,length=5), seq(0,100,length=5))
abline(0,1)
image(t(est3), zlim=c(-6,6), col=jet.colors(64), xaxt="n", yaxt="n", main="Estimate")
axis(1, seq(0,1,length=5), seq(0,100,length=5))
axis(2, seq(0,1,length=5), seq(0,100,length=5))
abline(0,1)




