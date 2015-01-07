# New pcox2 testing


load_all()
library(survival)
library(mgcv)
library(refundDevel)

data(sofa_fu)
data(sofa)

debug(pcox)


####################
# Parametric Terms #
####################

fit1 <- pcox(Surv(los, death) ~ age, data=sofa)
fit1a <- coxph(Surv(los, death) ~ age, data=sofa)




##################
# Smooth scalars #
##################

pdata <- data.frame(age=seq(min(sofa$age), max(sofa$age), by=.5))
fit2 <- pcox(Surv(los, death) ~ p(age, linear=FALSE, dbug=TRUE) + male, data=sofa)
fhat <- PredictMat(fit2$pcox$smooth[[1]][[1]], data=pdata) %*% fit2$coefficients[1:9]
qplot(pdata$age, fhat, geom="line")

# With simulated data
N <- 500
J <- 200
x <- runif(N, 0, 2*pi)
male <- rbinom(N, size = 1, prob=.5)
ftrue <- sin(x)
eta <- matrix(ftrue + .75*male, nrow=N, ncol=J)
Xdat <- data.frame(x=x, male=male)
data2 <- simTVSurv(eta, Xdat)
fit2.1 <- pcox(Surv(time, event) ~ p(x, linear=FALSE, dbug=TRUE) + male, data=data2)
pdata.1 <- data.frame(x=seq(0,2*pi,by=.1))
fhat.1 <- fit2.1$pcox$t.funcs[[1]](pdata.1) %*% fit2.1$coefficients[1:9]
qplot(pdata.1$x, fhat.1, geom="line") +
  geom_line(aes(y=sin(pdata.1$x)), col="red")



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
fit3.1 <- pcox(Surv(time,event) ~ bf(myX) + male, data=data3)
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



###################
# Historical TVCs #
###################
N <- 500
J <- 100
sind <- 1:J
X <- genX(N, seq(0,1,length=J))
beta <- makeBetaMat(J, genBeta1)

eta  <- sapply(1:J, function(j) {
  .75*male + X[,1:j,drop=F] %*% beta[j,1:j] / j
})
Xdat <- data.frame(myX=I(X), male=male)
data5 <- simTVSurv(eta, Xdat)
sind2 <- 1:ncol(data5$myX)

fit5.1 <- pcox(Surv(time,event) ~ male + hf(myX, sind = sind2, dbug=TRUE),
               data=data5)
est5.1 <- getHCEst(fit5.1$pcox$smooth[[1]][[1]], 1:J,
                   coefs = fit5.1$coefficients[-1])



