# New pcox2 testing

library(devtools)
dev_mode()
load_all()
library(survival)
library(mgcv)
library(refundDevel)

data(sofa_fu)
data(sofa)



# Take these out later?
library(pryr)
library(reshape2)

###########################
# Simple Parametric Model #
###########################

fit1 <- pcox(Surv(los, death) ~ age, data=sofa)
fit1a <- coxph(Surv(los, death) ~ age, data=sofa)
pre1.1 <- predict(fit1)
pre1.2 <- predict(fit1, newdata=sofa)
range(pre1.1 - pre1.2, na.rm=T) # Should be 0


##################
# Smooth scalars #
##################

pdata <- data.frame(age=seq(min(sofa$age), max(sofa$age), by=.5))
fit2 <- pcox(Surv(los, death) ~ p(age, linear=FALSE, dbug=TRUE) + male, data=sofa)
fhat <- PredictMat(fit2$pcox$smooth[[1]][[1]], data=pdata) %*% fit2$coefficients[1:9]
qplot(pdata$age, fhat, geom="line")

# STANDARD MODEL
set.seed(1235)
N <- 500
J <- 200
x <- runif(N, 0, 2*pi)
male <- rbinom(N, size = 1, prob=.5)
ftrue <- sin(x)
eta <- matrix(ftrue + .75*male, nrow=N, ncol=J)
Xdat <- data.frame(x=x, male=male)
data2 <- simTVSurv(eta, Xdat)
fit2.1 <- pcox(Surv(time, event) ~ p(x, linear=FALSE, dbug=TRUE) +
                 male, data=data2)

# Predict
pre2.1.1 <- predict(fit2.1)
pre2.1.2 <- predict(fit2.1, newdata=data2)
pre2.1.3 <- predict(fit2.1, newdata=data2[-4]) # Throws an error
range(pre2.1.1 - pre2.1.2, na.rm=T) # Should be 0

# Coef
pdata.1 <- data.frame(x=seq(0,2*pi,by=.1))
fhat.1 <- fit2.1$pcox$t.funcs[[1]](pdata.1) %*% fit2.1$coefficients[1:9]
qplot(pdata.1$x, fhat.1, geom="line") +
  geom_line(aes(y=sin(pdata.1$x)), col="red")


# TIME-VARYING MODEL
fit2.2 <- pcox(Surv(time, event) ~ p(x, linear=FALSE, tv=T) +
                 male, data=data2)




######################
# Baseline Functions #
######################

N <- 500
J <- 101
sind <- seq(0,1,length=J)
X <- genX(N, sind)
beta <- sin(2*pi*sind)
eta <- matrix((X%*%beta/J + .75*male), nrow=N, ncol=200) 
data3 <- simTVSurv(eta)
data3$myX <- I(X)
data3$male <- male
fit3.1 <- pcox(Surv(time,event) ~ bf(myX, bs="ps", sind=sind) + male, data=data3)

# Predict
pre3.1.1 <- predict(fit3.1)
pre3.1.2 <- predict(fit3.1, newdata=data3)
pre3.1.3 <- predict(fit3.1, newdata=data3[-3]) # Should throw an error
range(pre3.1.1 - pre3.1.2, na.rm=T) # Should be 0



# TO DO:  Try some things with sindices next


sm <- fit3.1$pcox$smooth[[1]][[1]]

pdata.31 <- data.frame(smat=I(matrix(sind, nrow=N, ncol=J)),
                       LX=I(matrix(1, nrow=N, ncol=J)))
pmat.31 <- PredictMat(sm, data = pdata.31)
pmat.32 <- fit3.1$pcox$t.funcs[[1]](data.frame(myX=I(X)))




pdata.3 <- data.frame(smat=sind, LX=1)
bhat3.1v1 <- PredictMat(fit3.1$pcox$smooth[[1]][[1]], pdata.3) %*%
  fit3.1$coef[-11]
bhat3.1 <- fit3.1$pcox$t.funcs[[1]](pdata.3) %*% fit3.1$coef[-11]

ggplot(pdata.3, mapping=aes(x = smat)) +
  geom_line(aes(y=beta)) +
  geom_line(aes(y=bhat3.1v1), col="red")



###################
# Concurrent TVCs #
###################

N <- 500
J <- 200
sind <- 1:J
X <- genX(N, seq(0,1,length=J))
eta <- t(sapply(1:N, function(i) {
  1.5*X[i,] + .75*male[i]
}))
Xdat <- data.frame(myX=I(X), male=male)
data4 <- simTVSurv(eta, Xdat)
sind2 <- 1:ncol(data4$myX)
fit4.1 <- pcox(Surv(time,event) ~ male + cf(myX, sind = sind2, dbug=TRUE),
               data=data4)

# Predict
pre4.1.1 <- predict(fit4.1)
pre4.1.2 <- predict(fit4.1, newdata=data4, stimes=data4$time)
pre4.1.3 <- predict(fit4.1, newdata=data4[-4])
range(pre4.1.1 - pre4.1.2, na.rm=T) # Should be 0



fit4.2 <- pcox(Surv(time,event) ~ male + cf(myX, tv=TRUE, sind = sind2, dbug=TRUE),
               data=data4)



###################
# Historical TVCs #
###################
N <- 500
J <- 100
sind <- 1:J
X <- genX(N, seq(0,1,length=J))
beta <- makeBetaMat(J, genBeta1)

eta  <- sapply(1:J, function(j) {
  .75*male + X[,1:j,drop=F] %*% beta[j,1:j]
})
Xdat <- data.frame(myX=I(X), male=male)
data5 <- simTVSurv(eta, Xdat)
sind2 <- 1:ncol(data5$myX)
fit5.1 <- pcox(Surv(time,event) ~ male + hf(myX, sind = sind2, dbug=TRUE),
               data=data5)

# Predict
pre5.1.1 <- predict(fit5.1)
pre5.1.2 <- predict(fit5.1, newdata=data5, stimes=data5$time)
pre5.1.3 <- predict(fit5.1, newdata=data5[-4], stimes=data5$time) # Should throw an error
range(pre5.1.1 - pre5.1.2, na.rm=T) # Should be 0 - confirms correct predictions for training data



# Coefficient
est5.1 <- getHCEst(fit5.1$pcox$smooth[[1]][[1]], 1:J,
                   coefs = fit5.1$coefficients[-1])

par(mfrow=c(1,2))
image.plot(t(beta))
image.plot(t(est5.1))

# Domain-standardized
fit5.1b <- pcox(Surv(time,event) ~ male +
                  hf(myX, sind = sind2, s.transform = "s/t", dbug=TRUE),
                data=data5)
pre5.1b <- predict(fit5.1b)








# Nonlinear

time5.2 <- system.time(fit5.2 <- pcox(Surv(time,event) ~ male + hf(myX, sind = sind2, dbug=TRUE, linear = F),
               data=data5))

