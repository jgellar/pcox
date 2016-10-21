library(pcox)
source(system.file("tests/setup-testthat.R", package = "pcox"))

# ------------------------------------------------------------------------------
context("basic functionality for parametric effects")
# h(x) = .6*x + .75*male 

test_that("parametric model works", {
  set.seed(121212)
  eta <- matrix(.6*x + .75*male, nrow=N, ncol=J)
  dat <- simTVSurv(eta, data.frame(x=x, male=male))
  fit <- pcox(Surv(time, event) ~ x + male, data=dat)
  fit_coxph <- survival::coxph(Surv(time, event) ~ x + male, data=dat)
  expect_equivalent(coef(fit), coef(fit_coxph))
  expect_equivalent(predict(fit),  predict(fit_coxph))
  expect_equivalent(predict(fit),  predict(fit, newdata=dat))
})
test_that("parametric model works with factors", {
  set.seed(121212)
  eta <- matrix(.6*x + .75*male, nrow=N, ncol=J)
  dat <- simTVSurv(eta, data.frame(x=x, male=male))
  dat$fmale <- factor(dat$male)
  fit <- pcox(Surv(time, event) ~ x + male, data=dat)
  fit_factor <- pcox(Surv(time, event) ~ x + fmale, data=dat)
  expect_equivalent(coef(fit), coef(fit_factor))
  expect_equivalent(predict(fit),  predict(fit_factor))
})
test_that("parametric model works with NAs", {
  set.seed(121212)
  eta <- matrix(.6*x + .75*male, nrow=N, ncol=J)
  dat <- simTVSurv(eta, data.frame(x=x, male=male))
  dat$xna <- dat$x; dat$xna[1:5] <- NA
  dat$timena <- dat$time; dat$timena[1:5] <- NA
  dat$eventna <- dat$event; dat$eventna[1:5] <- NA
  fit <- pcox(Surv(time, event) ~ x + male, data=dat)
  fit_xna <- pcox(Surv(time, event) ~ xna + male, data=dat)
  fit_eventna <- pcox(Surv(time, eventna) ~ x + male, data=dat)
  fit_timena <- pcox(Surv(timena, event) ~ x + male, data=dat)
  expect_equal(unname(coef(fit)), unname(coef(fit_xna)), tolerance = 0.1)
  expect_equal(unname(coef(fit)), unname(coef(fit_eventna)), tolerance = 0.1)
  expect_equal(unname(coef(fit)), unname(coef(fit_timena)), tolerance = 0.1)
})

# ------------------------------------------------------------------------------
context("basic functionality for nonlinear timeconstant effects")
# h(x) = sin(x) + .75*male 

test_that("smooth term works with mgcv options", {
  set.seed(121212)
  eta <- matrix(sin(x), nrow=N, ncol=J)
  dat <- simTVSurv(eta, data.frame(x=x))
  
  fit_8 <- pcox(Surv(time, event) ~ p(x, linear=FALSE, k=8), data=dat)
  expect_gt(cor(coef(fit_8)$value, sin(coef(fit_8))$x), .95)
  expect_equal(predict(fit_8), predict(fit_8, newdata=dat))
  
  fit_25 <- pcox(Surv(time, event) ~ p(x, linear=FALSE, k=25), data=dat)
  expect_gt(cor(coef(fit_25)$value, sin(coef(fit_8))$value), .95)
  expect_equal(fit_25$pcox$smooth[[1]]$bs.dim, 25)
  
  fit_ps <- pcox(Surv(time, event) ~ p(x, bs = "ps", linear=FALSE), data=dat)
  expect_gt(cor(coef(fit_ps)$value, sin(coef(fit_8))$value), .95)
  expect_s3_class(fit_ps$pcox$smooth[[1]], "pspline.smooth")
  
  # failing:
  # pcox(Surv(time, event) ~ p(x, basisype  = "te", linear=FALSE), data=dat)
  # pcox(Surv(time, event) ~ p(x, linear=T, k=8), data=dat)
  # pcox(Surv(time, event) ~ p(x, basisype  = "te"), data=dat) 
})

test_that("smooth term works with NA", {  
  set.seed(121212)
  eta <- matrix(sin(x), nrow=N, ncol=J)
  dat <- simTVSurv(eta, data.frame(x=x))
  dat$xna <- x; dat$xna[1:5] <- NA
  fit <- pcox(Surv(time, event) ~ p(x, linear=FALSE), data=dat)
  fit_na <- pcox(Surv(time, event) ~ p(xna, linear=FALSE), data=dat)
  expect_gt(cor(coef(fit)$value, sin(coef(fit_na))$x), .99)
  
  # failing:
  # expect_equal(predict(fit_na)[1:5], 1*rep(NA, 5))  # !!
  
})  

# ------------------------------------------------------------------------------
context("basic functionality for timevarying effects of scalars")
# h(t) = z * sin(t)

test_that("linear time-varying effect of scalar works", {
  set.seed(121212)
  eta <- matrix(sin(2 * pi*(1:J)/J) %x% male, nrow=N, ncol=J)
  dat <- simTVSurv(eta, data.frame(z=z, male = male, fmale=factor(male)))
  fit <- pcox(Surv(time, event) ~ p(male, linear=TRUE, tv=T), data=dat)
  est <- coef(fit)
  est_manual <- mgcv::PredictMat(fit$pcox$smooth[[1]], 
    data=data.frame(t=est$t, male=1)) %*% fit$coefficients
  expect_gt(cor(est$value, sin(2 * pi* est$t / J)), 0.94)
  expect_equal(predict(fit), predict(fit, newdata=dat, stimes = dat$time))
  expect_equal(est$value, est_manual)
  # failing:
  #fit <- pcox(Surv(time, event) ~ p(fmale, linear=TRUE, tv=T), data=dat)
})  

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


summary(fit2.1)
est2.1a <- coef(fit2.1)
est2.1b <- drop(fit2.1$pcox$t.funcs[[1]](est2.1a["x"]) %*%
    fit2.1$coefficients[1:7])
ggplot(est2.1a, aes(x, value)) + geom_line(colour="red", size=2) +
  geom_line(aes(y=sin(x)), size=2) +
  geom_line(aes(y=est2.1b), colour="blue", linetype="dashed", size=2)
pre2.1a <- predict(fit2.1)
pre2.1b <- predict(fit2.1, newdata=dat2.1)
range(pre2.1a - pre2.1b, na.rm=T)


