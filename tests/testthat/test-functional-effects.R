library(pcox)

set.seed(12354)
N <- 500
J <- 200
x <- runif(N, 0, 2*pi)
z <- rnorm(N)
male <- rbinom(N, size = 1, prob=.5)
somefactor <- factor(sample(1:3, N, replace = TRUE))
sind <- seq(0,1,length=J)
X <- genX(N, sind)
K <- 100
Z <- genX(N, seq(0,1,length=K))
L <- 100

# - check NA handling (#11)
# - check predict, coef
# 

# ------------------------------------------------------------------------------
context("estimate timeconstant linear baseline functional effects")
# beta <- sin(2*pi*sind)
# eta3.1 <- matrix((X%*%sin(2*pi*sind)/J + .75*male), nrow=N, ncol=200)
# dat3.1 <- simTVSurv(eta3.1)
# dat3.1$myX <- X
# dat3.1$male <- male
# fit3.1 <- pcox(Surv(time,event) ~ bf(myX, bs="ps", sind=sind) + male,
#   data=dat3.1)
# est3.1a <- coef(fit3.1)
# est3.1b <- mgcv::PredictMat(fit3.1$pcox$smooth[[1]],
#   data=data.frame(myX.smat=est3.1a$s, myX.LX=1)) %*%
#   fit3.1$coefficients[1:10]
# ggplot(est3.1a, aes(s, value)) + geom_line(colour="red", size=2) +
#   geom_line(aes(y=sin(2*pi*s)), size=2) +
#   geom_line(aes(y=est3.1b), colour="blue", linetype="dashed", size=2)
# pre3.1a <- predict(fit3.1)
# pre3.1b <- predict(fit3.1, newdata=dat3.1)
# range(pre3.1a - pre3.1b, na.rm=T) # Should be 0
# summary(fit3.1)