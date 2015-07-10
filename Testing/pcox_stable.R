# Stable pcox tests: Run this first!


library(devtools)
dev_mode()
load_all()
#library(survival)
#library(mgcv)

library(ggplot2)
library(gridExtra)
library(RColorBrewer)

# Take these out later?
#library(pryr)
library(reshape2)



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
L <- 100

plotMe <- function(est, lims=range(est$value)) {
  ggplot(est, aes(s, t)) + 
    geom_tile(aes(fill=value, colour=value)) +
    theme_bw() +
    scale_fill_gradientn(name="", limits=lims,
                         colours=rev(brewer.pal(11,"Spectral"))) +
    scale_colour_gradientn(name="", limits=lims,
                           colours=rev(brewer.pal(11,"Spectral"))) +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0))
}


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
fit1.1

time1 <- dat1.1$time
evnt1 <- dat1.1$event
fit1.2 <- pcox(Surv(time1, evnt1) ~ x + male, data=dat1.1)
fit1.2
summary(fit1.2)


##################
# Smooth scalars #
##################

# STANDARD MODEL, i.e. ~ \beta(x)
eta2.1 <- matrix(sin(x) + .75*male, nrow=N, ncol=J)
dat2.1 <- simTVSurv(eta2.1, data.frame(x=x, male=male))
fit2.1 <- pcox(Surv(time, event) ~ p(x, linear=FALSE, k=8) +
                 male, data=dat2.1)
fit2.1
summary(fit2.1)
est2.1a <- coef(fit2.1)
est2.1b <- drop(fit2.1$pcox$t.funcs[[1]](est2.1a["x"]) %*%
                  fit2.1$coefficients[1:7])
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
est2.2b <- mgcv::PredictMat(fit2.2$pcox$smooth[[1]],
                      data=data.frame(t=est2.2a$t, x=1)) %*% 
  fit2.2$coefficients[1:10]
ggplot(est2.2a, aes(t, value)) + geom_line(colour="red", size=2) +
  geom_line(aes(y=sin(2*pi*t/J)), size=2) +
  geom_line(aes(y=est2.2b), colour="blue", linetype="dashed", size=2)
pre2.2a <- predict(fit2.2)
pre2.2b <- predict(fit2.2, newdata=dat2.2, stimes=dat2.2$time)
range(pre2.2a - pre2.2b, na.rm=T) # Should be 0
fit2.2
summary(fit2.2)

# Two time-varying terms
eta2.3 <- matrix(sin(2*pi*(1:J)/J) %x% z, nrow=N, ncol=J) +
  matrix(cos(2*pi*(1:J)/J) %x% male, nrow=N, ncol=J)
dat2.3 <- simTVSurv(eta2.3, data.frame(x=z, male=male))
fit2.3 <- pcox(Surv(time, event) ~ p(x, linear=TRUE, tv=T) +
                 p(male, linear=T, tv=T), data=dat2.3)
fit2.3
summary(fit2.3)
est2.3a <- coef(fit2.3)
est2.3b <- coef(fit2.3, select=2)
p1 <- ggplot(est2.3a, aes(t, value)) + geom_line(colour="red", size=2) +
  geom_line(aes(y=sin(2*pi*t/J)), size=2) + ylim(c(-1.5,1.5))
p2 <- ggplot(est2.3b, aes(t, value)) + geom_line(colour="red", size=2) +
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
est3.1b <- mgcv::PredictMat(fit3.1$pcox$smooth[[1]],
                      data=data.frame(myX.smat=est3.1a$s, myX.LX=1)) %*% 
  fit3.1$coefficients[1:10]
ggplot(est3.1a, aes(s, value)) + geom_line(colour="red", size=2) +
  geom_line(aes(y=sin(2*pi*s)), size=2) +
  geom_line(aes(y=est3.1b), colour="blue", linetype="dashed", size=2)
pre3.1a <- predict(fit3.1)
pre3.1b <- predict(fit3.1, newdata=dat3.1)
range(pre3.1a - pre3.1b, na.rm=T) # Should be 0

# TIME-VARYING BASELINE FUNCTION
# \beta(s,t) = 5*sin(2*pi*s + pi*t/100)
# \eta_i(t) = \int X_i(s)\beta(s,t)ds + .75*male_i

beta3.2 <- expand.grid(s=sind, t=(0:L))
beta3.2$value <- 5*sin(2*pi*beta3.2$s + pi*beta3.2$t/100)
#*cos(2*pi*beta3.2$t/100)
beta3.2mat <- acast(beta3.2, t ~ s)
eta3.2 <- t(apply(X, 1, function(x.i) {
  x.i %*% t(beta3.2mat)/J
})) + .75*matrix(male, nrow=N, ncol=(L+1))
dat3.2 <- simTVSurv(eta3.2)
dat3.2$myX <- X
dat3.2$male <- male
fit3.2 <- pcox(Surv(time,event) ~ bf(myX, bs="ps", basistype="te", sind=sind,
                                     tv=TRUE) + male, data=dat3.2)
est3.2 <- coef(fit3.2)
p3.2a <- plotMe(beta3.2, c(-6,6))
p3.2b <- plotMe(est3.2,  c(-6,6))
p3.2c <- arrangeGrob(p3.2a, p3.2b, nrow=1)
p3.2c

# ADDITIVE BASELINE (here I'm using t instead of x just so plotMe works)
beta3.3 <- expand.grid(t=seq(min(X), max(X), length=J), s=sind)
beta3.3$value <- 20*sin(2*pi*beta3.3$s + pi*beta3.3$t/20)
beta3.3mat <- acast(beta3.3, t ~ s)
eta3.3 <- apply(X, 1, function(x.i) {
  mean(20*sin(2*pi*sind + pi*x.i/20))
}) + .75*male
dat3.3 <- simTVSurv(matrix(eta3.3, nrow=N, ncol=L))
dat3.3$myX <- X
dat3.3$male <- male
fit3.3 <- pcox(Surv(time,event) ~ bf(myX, bs="ps", basistype="te", sind=sind,
                                     linear=FALSE) + male, data=dat3.3)
est3.3 <- coef(fit3.3)
names(est3.3)[2] <- "t"
p3.3a <- plotMe(beta3.3, c(-30,30))
p3.3b <- plotMe(est3.3,  c(-30,30))
p3.3c <- arrangeGrob(p3.3a, p3.3b, nrow=1)
p3.3c





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
fit4.1
summary(fit4.1)


# Time-varying version
eta4.2 <- t(sapply(1:N, function(i) { sin(2*pi*(1:J)/J)*X[i,] + .75*male[i] }))
dat4.2 <- simTVSurv(eta4.2, data.frame(myX=I(X), male=male))
fit4.2 <- pcox(Surv(time,event) ~ male + cf(myX, tv=TRUE,
                                            sind = (1:ncol(dat4.2$myX))),
               data=dat4.2)
fit4.2
summary(fit4.2)
est4.2a <- coef(fit4.2)
est4.2b <- mgcv::PredictMat(fit4.2$pcox$smooth[[1]],
                      data=data.frame(t=est4.2a$t, myX.t=1)) %*% 
  fit4.2$coefficients[2:11]
ggplot(est4.2a, aes(t, value)) + geom_line(colour="red", size=2) +
  geom_line(aes(y=sin(2*pi*t/J)), size=2) +
  geom_line(aes(y=est4.2b), colour="blue", linetype="dashed", size=2)
pre4.2a <- predict(fit4.2)
pre4.2b <- predict(fit4.2, newdata = dat4.2, stimes=dat4.2$time)


# Lagged concurrent TVC's
lag=5
eta4.3 <- 1.5*cbind(matrix(0, nrow=N, ncol=lag), X[, 1:(J-lag)]) +
  .75*matrix(male, nrow=N, ncol=J)
dat4.3 <- simTVSurv(eta4.3, data.frame(myX=I(X), male=male))
fit4.3 <- pcox(Surv(time,event) ~ male + cf(myX, sind = (1:ncol(dat4.3$myX)), lag=lag),
               data=dat4.3)
fit4.3
summary(fit4.3)
est4.3 <- coef(fit4.3)
pre4.3a <- predict(fit4.3)
pre4.3b <- predict(fit4.3, newdata=dat4.3, stimes=dat4.3$time)
range(pre4.3a - pre4.3b, na.rm=T) # Should be 0


###################
# Historical TVCs #
###################


beta <- makeBetaMat(K, genBeta1)
eta5.1  <- sapply(1:K, function(k) {
  .75*male + (Z[,1:k,drop=F] %*% beta[k,1:k])/k
})
dat5.1 <- simTVSurv(eta5.1, data.frame(myX=I(Z), male=male))
sinds <- 1:ncol(dat5.1$myX)

fit5.1 <- pcox(Surv(time,event) ~ male + hf(myX, sind = sinds), data=dat5.1)
fit5.2 <- pcox(Surv(time,event) ~ male + hf(myX, sind = sinds, basistype="te",
                                            bs="ps"), data=dat5.1)
fit5.3 <- pcox(Surv(time,event) ~ male + hf(myX, sind = sinds,
                                            transform="lagged"), data=dat5.1)
fit5.4 <- pcox(Surv(time,event) ~ male + hf(myX, sind = sinds,
                                            transform="standardized"), data=dat5.1)
fit5.5 <- pcox(Surv(time,event) ~ male + hf(myX, sind = sinds, transform="standardized",
                                            basistype = "te", bs="ps"), data=dat5.1)
fit5.6 <- pcox(Surv(time,event) ~ male + hf(myX, sind = sinds, transform="linear",
                                            bs="ps"), data=dat5.1)
fit5.7 <- pcox(Surv(time,event) ~ male + hf(myX, sind = sinds, transform="quadratic",
                                            bs="ps"), data=dat5.1)
fit5.8 <- pcox(Surv(time,event) ~ male + hf(myX, sind = sinds, transform="noInteraction",
                                            bs="ps"), data=dat5.1)
fit5.1
summary(fit5.1)
plotMe(coef(fit5.1), c(-6,6))
plotMe(coef(fit5.2), c(-6,6))
plotMe(coef(fit5.3), c(-6,6))
plotMe(coef(fit5.4), c(-6,6))
plotMe(coef(fit5.5), c(-6,6))
plotMe(coef(fit5.6), c(-6,6))
plotMe(coef(fit5.7), c(-6,6))
plotMe(coef(fit5.8), c(-6,6))

est5.1 <- coef(fit5.1)
truebeta <- genBeta1(est5.1$s, est5.1$t)
amse <- sapply(1:8, function(i) {
  est <- coef(get(paste0("fit5.",i)))
  mean((truebeta - est$value)^2)
})
amse
plot(amse)



# Special Terms from coxph

# Strata and Cluster
bladder1 <- bladder[bladder$enum < 5, ] 
fit.S1 <- coxph(Surv(stop, event) ~ (rx + size + number) * strata(enum) + 
                  cluster(id), bladder1)
fit.S2 <- pcox( Surv(stop, event) ~ (rx + size + number) * strata(enum) + 
                  cluster(id), bladder1)
all.equal(coef(fit.S1), coef(fit.S2))

# tt
fit.tt1 <- coxph(Surv(time, status) ~ ph.ecog + tt(age), data=lung,
                 tt=function(x,t,...) pspline(x + t/365.25))
fit.tt2 <- pcox( Surv(time, status) ~ ph.ecog + tt(age), data=lung,
                 tt=function(x,t,...) pspline(x + t/365.25))
all.equal(coef(fit.tt1), coef(fit.tt2))

# Multiple tt's (some from user, some from p())
fit.tt3 <- pcox( Surv(time, status) ~ p(meal.cal, tv = TRUE) + sex + tt(age),
                 data=lung, tt=function(x,t,...) pspline(x + t/365.25))
fit.tt4 <- pcox( Surv(time, status) ~ tt(age) + p(meal.cal, tv = TRUE) + sex,
                 data=lung, tt=function(x,t,...) pspline(x + t/365.25))
all.equal(fit.tt3$coefficients,
          fit.tt4$coefficients[names(fit.tt3$coefficients)])

# coxph.penalty terms
fit.ridge1 <- coxph(Surv(futime, fustat) ~ rx + ridge(age, ecog.ps, theta=1),
                    ovarian)
fit.ridge2 <- pcox( Surv(futime, fustat) ~ rx + ridge(age, ecog.ps, theta=1),
                    ovarian)
all.equal(fit.ridge1$coefficients, fit.ridge2$coefficients)

fit.ps1 <- coxph(Surv(time, status) ~ ph.ecog + pspline(age,8), cancer)
fit.ps2 <- pcox( Surv(time, status) ~ ph.ecog + pspline(age,8), cancer)
all.equal(fit.ps1$coefficients, fit.ps2$coefficients)

fit.frailty1 <- coxph(Surv(time, status) ~ age + frailty(inst, df=4), lung)
fit.frailty2 <- pcox( Surv(time, status) ~ age + frailty(inst, df=4), lung)
all.equal(fit.frailty1$coefficients, fit.frailty2$coefficients)

