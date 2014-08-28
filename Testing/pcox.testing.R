# pcox.testing.R

# s(x)

N <- 500
J <- 101
x <- rnorm(N)
eta <- sin(x)
data1 <- simTVSurv(matrix(eta, nrow=N, ncol=J), Xdat=data.frame(x=x))

fit <- pcox(Surv(data1$time, data1$event) ~ s(data1$x))
fit <- pcox(Surv(time, event) ~ s(x), data=data1)


