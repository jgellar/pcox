# Stable pcox tests: Run this first!


library(devtools)
dev_mode()
load_all()
library(survival)
library(mgcv)
library(ggplot2)
library(gridExtra)


# Take these out later?
library(pryr)
library(reshape2)


library(ggplot2)
library(RColorBrewer)

# Set up variables

set.seed(12354)
N <- 500
J <- 200
x <- runif(N, 0, 2*pi)
z <- rnorm(N)
male <- rbinom(N, size = 1, prob=.5)
sind <- seq(0,1,length=J)
X <- genX(N, sind)
K <- 100
Z <- genX(N, seq(0,1,length=K))


###########################
# Simple Parametric Model #
###########################


eta1.1 <- matrix(.6*x + .75*male, nrow=N, ncol=J)
dat1.1 <- simTVSurv(eta1.1, data.frame(x=x, male=male))
fit1.1 <- pcox(Surv(time, event) ~ x + male, data=dat1.1)
est1.1 <- coef(fit1.1)
pre1.1a <- predict(fit1.1)
pre1.1b <- predict(fit1.1, newdata=dat1.1)
range(pre1.1a - pre1.1b)


##################
# Smooth scalars #
##################

# STANDARD MODEL, i.e. ~ \beta(x)
eta2.1 <- matrix(sin(x) + .75*male, nrow=N, ncol=J)
dat2.1 <- simTVSurv(eta2.1, data.frame(x=x, male=male))
fit2.1 <- pcox(Surv(time, event) ~ p(x, linear=FALSE, dbug=TRUE) +
                 male, data=dat2.1)
est2.1a <- coef(fit2.1)
est2.1b <- drop(fit2.1$pcox$t.funcs[[1]](est2.1a["x"]) %*%
                  fit2.1$coefficients[1:9])
ggplot(est2.1a, aes(x, value)) + geom_line(colour="red", size=2) +
  geom_line(aes(y=sin(x)), size=2) +
  geom_line(aes(y=est2.1b), colour="blue", linetype="dashed", size=2)
pre2.1a <- predict(fit2.1)
pre2.1b <- predict(fit2.1, newdata=dat2.1)
range(pre2.1a - pre2.1b, na.rm=T) # Should be 0

# LINEAR, TIME-VARYING MODEL, i.e. ~ \beta(t)*z
eta2.2 <- matrix(sin(2*pi*(1:J)/J) %x% z, nrow=N, ncol=J) +
  matrix(0.75*male, nrow=N, ncol=J)
dat2.2 <- simTVSurv(eta2.2, data.frame(x=z, male=male))
fit2.2 <- pcox(Surv(time, event) ~ p(x, linear=TRUE, tv=T) +
                 male, data=dat2.2)
est2.2a <- coef(fit2.2)
est2.2b <- PredictMat(fit2.2$pcox$smooth[[1]],
                      data=data.frame(t=est2.2a$t, x=1)) %*% 
  fit2.2$coefficients[1:10]
ggplot(est2.2a, aes(t, value)) + geom_line(colour="red", size=2) +
  geom_line(aes(y=sin(2*pi*t/J)), size=2) +
  geom_line(aes(y=est2.2b), colour="blue", linetype="dashed", size=2)
pre2.2a <- predict(fit2.2)
pre2.2b <- predict(fit2.2, newdata=dat2.2, stimes=dat2.2$time)
range(pre2.2a - pre2.2b, na.rm=T) # Should be 0

# Two time-varying terms
eta2.3 <- matrix(sin(2*pi*(1:J)/J) %x% z, nrow=N, ncol=J) +
  matrix(cos(2*pi*(1:J)/J) %x% male, nrow=N, ncol=J)
dat2.3 <- simTVSurv(eta2.3, data.frame(x=z, male=male))
fit2.3 <- pcox(Surv(time, event) ~ p(x, linear=TRUE, tv=T) +
                 p(male, linear=T, tv=T), data=dat2.3)
est2.3a <- coef(fit2.3)
p1 <- ggplot(est2.3a[[1]], aes(t, value)) + geom_line(colour="red", size=2) +
  geom_line(aes(y=sin(2*pi*t/J)), size=2) + ylim(c(-1.5,1.5))
p2 <- ggplot(est2.3a[[2]], aes(t, value)) + geom_line(colour="red", size=2) +
  geom_line(aes(y=cos(2*pi*t/J)), size=2) + ylim(c(-1.5,1.5))
p  <- arrangeGrob(p1,p2,nrow=1)
p



##################################
# Baseline Functional Predictors #
##################################

beta <- sin(2*pi*sind)
eta3.1 <- matrix((X%*%sin(2*pi*sind)/J + .75*male), nrow=N, ncol=200) 
dat3.1 <- simTVSurv(eta3.1)
dat3.1$myX <- X
dat3.1$male <- male
fit3.1 <- pcox(Surv(time,event) ~ bf(myX, bs="ps", sind=sind) + male,
               data=dat3.1)
est3.1a <- coef(fit3.1)
est3.1b <- PredictMat(fit3.1$pcox$smooth[[1]],
                      data=data.frame(myX.smat=est3.1a$s, myX.LX=1)) %*% 
  fit3.1$coefficients[1:10]
ggplot(est3.1a, aes(s, value)) + geom_line(colour="red", size=2) +
  geom_line(aes(y=sin(2*pi*s)), size=2) +
  geom_line(aes(y=est3.1b), colour="blue", linetype="dashed", size=2)
pre3.1a <- predict(fit3.1)
pre3.1b <- predict(fit3.1, newdata=dat3.1)
range(pre3.1a - pre3.1b, na.rm=T) # Should be 0


###################
# Concurrent TVCs #
###################

eta4.1 <- t(sapply(1:N, function(i) { 1.5*X[i,] + .75*male[i] }))
dat4.1 <- simTVSurv(eta4.1, data.frame(myX=I(X), male=male))
fit4.1 <- pcox(Surv(time,event) ~ male + cf(myX, sind = (1:ncol(dat4.1$myX))),
               data=dat4.1)
est4.1 <- coef(fit4.1)
pre4.1a <- predict(fit4.1)
pre4.1b <- predict(fit4.1, newdata=dat4.1, stimes=dat4.1$time)
range(pre4.1a - pre4.1b, na.rm=T) # Should be 0


# Time-varying version
eta4.2 <- t(sapply(1:N, function(i) { sin(2*pi*(1:J)/J)*X[i,] + .75*male[i] }))
dat4.2 <- simTVSurv(eta4.2, data.frame(myX=I(X), male=male))
fit4.2 <- pcox(Surv(time,event) ~ male + cf(myX, tv=TRUE,
                                            sind = (1:ncol(dat4.2$myX))),
               data=dat4.2)
est4.2a <- coef(fit4.2)
est4.2b <- PredictMat(fit4.2$pcox$smooth[[1]],
                      data=data.frame(t=est4.2a$t, myX.t=1)) %*% 
  fit4.2$coefficients[2:11]
ggplot(est4.2a, aes(t, value)) + geom_line(colour="red", size=2) +
  geom_line(aes(y=sin(2*pi*t/J)), size=2) +
  geom_line(aes(y=est4.2b), colour="blue", linetype="dashed", size=2)
pre4.2a <- predict(fit4.2)
pre4.2b <- predict(fit4.2, newdata = dat4.2, stimes=dat4.2$time)


###################
# Historical TVCs #
###################
beta <- makeBetaMat(K, genBeta1)

eta5.1  <- sapply(1:K, function(k) {
  .75*male + (Z[,1:k,drop=F] %*% beta[k,1:k])/k
})
dat5.1 <- simTVSurv(eta5.1, data.frame(myX=I(Z), male=male))
fit5.1 <- pcox(Surv(time,event) ~ male + hf(myX, sind = (1:ncol(dat5.1$myX))),
               data=dat5.1)
pre5.1a <- predict(fit5.1)
pre5.1b <- predict(fit5.1, newdata=dat5.1, stimes=dat5.1$time)
range(pre5.1a - pre5.1b, na.rm=T) # Should be 0 - confirms correct predictions for training data
est5.1 <- coef(fit5.1)
est5.1_old <- getHCEst(fit5.1$pcox$smooth[[1]], 1:K,
                       coefs = fit5.1$coefficients[-1])
lims <- range(est5.1$value)
ggplot(est5.1, aes(s, t)) + 
  geom_tile(aes(fill=value, colour=value)) +
  theme_bw() +
  scale_fill_gradientn(name="", limits=lims,
                       colours=rev(brewer.pal(11,"Spectral"))) +
  scale_colour_gradientn(name="", limits=lims,
                         colours=rev(brewer.pal(11,"Spectral"))) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0))

library(fields)
par(mfrow=c(1,2))
image.plot(t(beta), zlim=c(-6,6))
image.plot(t(est5.1_old), zlim=c(-6,6))

