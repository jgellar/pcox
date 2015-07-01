library(pcox)
library(pryr)
library(reshape2)
library(testthat)

context("Smooth scalar fits")


set.seed(1235)
N <- 500
J <- 200
x <- runif(N, 0, 2*pi)
male <- rbinom(N, size = 1, prob=.5)
ftrue <- sin(x)
eta <- matrix(ftrue + .75*male, nrow=N, ncol=J)
Xdat <- data.frame(x=x, male=male)
data2 <- simTVSurv(eta, Xdat)
fit2.1 <- pcox(Surv(time, event) ~ p(x, linear=FALSE) +
                 male, data=data2)

test_that("Base model fits", {
  expect_is(fit2.1, "pcox")
  expect_is(fit2.1, "coxph")
})


# test_that("coef method works", {
#   est2.1 <- coef(fit2.1)
#   expect_equal(dim(est2.1), c(101,4))
#   expect_equal(names(est2.1), c("value", "se", "se2", "x"))
#   fhat <- drop(fit2.1$pcox$t.funcs[[1]](est2.1["x"]) %*%
#                  fit2.1$coefficients[1:9])
#   expect_equal(est2.1$value, fhat)
#   eps <- est2.1$value - sin(est2.1$x)
#   expect_true(abs(mean(eps)) < .1)
#   expect_true(mean(abs(eps)) < .25)
#   expect_true(max(abs(eps)) < .5)
# })


test_that("predict method works", {
  pre2.1.1 <- predict(fit2.1)
  pre2.1.2 <- predict(fit2.1, newdata=data2)
  expect_equal(pre2.1.1, pre2.1.2)
  expect_error(predict(fit2.1, newdata=data2[-4]))
})


# test_that("Time-varying model fits and methods still work", {
#   x <- rnorm(N)
#   ftrue <- sin(2*pi*(1:J)/J)
#   eta <- matrix(ftrue %x% x, nrow=N, ncol=J) + matrix(0.75*male, nrow=N, ncol=J)
#   Xdat <- data.frame(x=x, male=male)
#   data2 <- simTVSurv(eta, Xdat)
#   fit2.2 <- pcox(Surv(time, event) ~ p(x, linear=TRUE, tv=T) +
#                    male, data=data2)
#   expect_is(fit2.2, "pcox")
#   
#   # coef Method
#   est2.2 <- coef(fit2.2)
#   fhat <- drop(fit2.2$pcox$t.funcs[[1]](rep(1,nrow(est2.2)), est2.2[["t"]]) %*%
#                  fit2.2$coefficients[1:10])
#   expect_equal(est2.2$value, fhat)
#   eps <- est2.2$value - sin(est2.2$t*2*pi/J)
#   expect_true(abs(mean(eps)) < 0.1)
#   expect_true(mean(abs(eps)) < 0.15)
#   
#   # predict method
#   pre2.2 <- predict(fit2.2)
#   
#   
#   pre2.1.1 <- predict(fit2.1)
#   pre2.1.2 <- predict(fit2.1, newdata=data2)
#   #expect_equal(pre2.1.1, pre2.1.2)
#   expect_error(predict(fit2.1, newdata=data2[-4]))
#   
#   
#   
#   
# })


