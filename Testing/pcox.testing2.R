# Some of these tests require running pcox_stable.R first.

data(sofa)


# Simple parameteric model

fit1 <- pcox(Surv(los, death) ~ age, data=sofa)
fit1a <- coxph(Surv(los, death) ~ age, data=sofa)
pre1.1 <- predict(fit1)
pre1.2 <- predict(fit1, newdata=sofa)
range(pre1.1 - pre1.2, na.rm=T) # Should be 0


fit2 <- pcox(Surv(los, death) ~ age*male + Charlson, data=sofa)
fit2a <- coxph(Surv(los, death) ~ age*male + Charlson, data=sofa)



##################
# Smooth scalars #
##################


# Random effect (frailty)
fit2.4 <- pcox(Surv(time, event) ~ p(x, linear=FALSE, k=8) + male +
                frailty(rep(1:4, N/4)),
               data=dat2.1)



fit2 <- pcox(Surv(los, death) ~ p(age, linear=FALSE) + male, data=sofa)
pdata <- data.frame(age=seq(min(sofa$age), max(sofa$age), by=.5))
fhat <- PredictMat(fit2$pcox$smooth[[1]], data=pdata) %*% fit2$coefficients[1:9]
qplot(pdata$age, fhat, geom="line")



pre2.1c <- predict(fit2.1, newdata=dat2.1[-3]) # Should throw an error
pre2.2c <- predict(fit2.2, newdata=dat2.2[-3],
                   stimes=dat2.2$time) # Should throw an error

# VARIABLE-COEFFICIENT SCALAR TERMS
tmp <- pcox(Surv(time,event) ~ p(male, by=x, linear=TRUE, tv=FALSE),
            data=data2)




##################################
# Baseline Functional Predictors #
##################################

pre3.1.3 <- predict(fit3.1, newdata=data3[-3]) # Should throw an error

# TO DO: Try some things with sindices next
# TO DO: Variable-domain functional predictors


###################
# Concurrent TVCs #
###################

pre4.1c <- predict(fit4.1, newdata=dat4.1[-4])
pre4.2c <- predict(fit4.2, newdata = dat4.2[-4], stimes=dat4.2$time)







###################
# Historical TVCs #
###################


pre5.1c <- predict(fit5.1, newdata=dat5.1[-4],
                   stimes=dat5.1$time) # Should throw an error

# Doesn't work yet:
ggplot(est5.1, aes(myX.s, myX.t)) + 
  geom_triangle(value)
par(mfrow=c(1,2))
image.plot(t(beta), zlim=c(-6,6))
image.plot(t(est5.1_old), zlim=c(-6,6))

# Domain-standardized
fit5.1a <- pcox(Surv(time,event) ~ male +
                hf(myX, sind = sind2, s.transform = "s/t", dbug=TRUE),
                data=data5)
pre5.1a <- predict(fit5.1a)

# Lagged model
fit5.1b <- pcox(Surv(time,event) ~ male +
                hf(myX, sind = sind2, s.transform = "s-t", dbug=TRUE),
                data=data5)
pre5.1b <- predict(fit5.1b)

# Linear interaction
fit5.1c <- pcox(Surv(time,event) ~ male +
                hf(myX, sind = sind2, s.transform = "s/t", dbug=TRUE,
                   bs="pb", xt="linear"),
                data=data5)
pre5.1c <- predict(fit5.1c)
# Nonlinear
time5.2 <- system.time(fit5.2 <- pcox(Surv(time,event) ~ male + hf(myX, sind = sind2, dbug=TRUE, linear = F),
                                      data=data5))



# Additive, Time-varying Baseline Function (TRIVARIATE!)
# F(x,s,t) = sin(pi*x/20 + 2pi*s - pi*t/100)
beta3.4 <- expand.grid(x=seq(min(X), max(X), length=50), s=sind,
                       t=(0:L))
beta3.4$value <- 20*sin(2*pi*beta3.4$s + pi*beta3.4$x/20 - pi*beta3.4$t/100)
beta3.4mat <- acast(beta3.4, x ~ t ~ s\  
eta3.4 <- sapply(0:L, function(l) {
  apply(X, 1, function(x.i) {
    mean(20*sin(2*pi*sind + pi*x.i/20 + pi*l/100))
  }) + .75*male
})
dat3.4 <- simTVSurv(matrix(eta3.4, nrow=N, ncol=L))
dat3.4$myX <- X
dat3.4$male <- male
fit3.4 <- pcox(Surv(time,event) ~ bf(myX, sind=sind,
                                     linear=FALSE, tv=TRUE) + male,
               data=dat3.4)
est3.4 <- coef(fit3.4)
est3.4 <- mycoef(fit3.4)

names(est3.4)[2] <- "t"
p3.4a <- plotMe(beta3.4, c(-30,30))
p3.4b <- plotMe(est3.4,  c(-30,30))
p3.4c <- arrangeGrob(p3.4a, p3.4b, nrow=1)
p3.4c



