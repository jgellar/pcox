# New pcox2 testing


load_all()
library(survival)
library(mgcv)
library(refundDevel)

data(sofa_fu)
data(sofa)


fit1 <- pcox2(Surv(time,event) ~ age, data=sofa_fu)

fit2 <- pcox2(Surv(time,event) ~ s(age, additive=TRUE, dbug=TRUE), data=sofa_fu)


N <- 500
J <- 200
x1 <- rbinom(N, size = 1, .5)
x2 <- runif(N, 0, 2*pi)
Xdat <- data.frame(x1=x1, x2=x2)
xb <- 3*sin(x2)
eta <- matrix(xb, nrow=N, ncol=J)

data2 <- simTVSurv(eta, Xdat)
fit.1 <- survreg(Surv(time,event) ~ pspline(x2), data=data2)
fh1 <- predict(fit.1, type="terms")[,1]

fit.2 <- pcox2(Surv(time,event) ~ s(x2, additive=TRUE, dbug=TRUE), data=data2)
fh2 <- fit.2$pcox$tt[[1]](data2$x2,data2$x2) %*% fit.2$coef

x.seq <- seq(0,2*pi,by=.1)
plot(x.seq, 3*sin(x.seq), type="l", lwd=3, ylim=c(-4,4))
points(data2$x2, fh1, col="blue")
points(data2$x2, fh2, col="red")



xb <- 2*x1 + 3*sin(x2)
eta <- matrix(xb, nrow=N, ncol=J)
data3 <- simTVSurv(eta, Xdat)
fit.3 <- pcox(Surv(time,event) ~ x1 + sm(x2, additive=TRUE, dbug=TRUE), data=data3)
fh3 <- fit.3$pcox$tt[[2]](data2$x2,data2$x2) %*% fit.3$coef[-1]

points(Xdat$x2, fh3, col="purple")


head(fit.3$coef)



xb <- 2*x1 + 3*sin(x2)
eta <- matrix(xb, nrow=N, ncol=J)
data3 <- simTVSurv(eta, Xdat)
fit.4 <- pcox2(Surv(time,event) ~ x1 + sm(x2, additive=TRUE, tv=TRUE, dbug=TRUE, k=20),
               data=data3)


newx <- seq(0, 2*pi, by=.1)
newt <- min(data3$time):max(data3$time)
newxmat <- matrix(newx, nrow=length(newx), ncol=length(newt))
newtmat <- matrix(newt, nrow=length(newx), ncol=length(newt), byrow=T)

fh4 <- fit.4$pcox$tt[[2]](as.vector(newxmat),as.vector(newtmat)) %*% fit.4$coef[-1]
fh4.mat <- matrix(fh4, nrow=length(newx), ncol=length(newt))


fit3 <- pcox(Surv(time, event) ~ )


