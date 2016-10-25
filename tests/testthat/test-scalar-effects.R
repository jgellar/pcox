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


# ------------------------------------------------------------------------------
context("estimate scalar parametric effects")
# h(x) = .6*x + .75*male 

test_that("parametric model works", {
  set.seed(121212)
  eta <- matrix(.6*x + .75*male, nrow=N, ncol=J)
  dat <- simTVSurv(eta, data.frame(x=x, male=male))
  fit <- pcox(Surv(time, event) ~ x + male, data=dat)
  fit_coxph <- survival::coxph(Surv(time, event) ~ x + male, data=dat)
  expect_equivalent(coef(fit), coef(fit_coxph))
  expect_equivalent(predict(fit), predict(fit_coxph))
  expect_equivalent(predict(fit), predict(fit, newdata=dat))
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
context("estimate nonlinear timeconstant scalar effects")
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
  fit_na <- pcox(Surv(time, event) ~ p(xna, linear=FALSE), data=dat)
  expect_gt(cor(coef(fit_na)$value, sin(coef(fit_na))$x), .99)
  
  # failing:
  # expect_equal(predict(fit_na)[1:5], 1*rep(NA, 5))  # !!
})  

# ------------------------------------------------------------------------------
context("estimate timevarying scalar effects")
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
  expect_error(pcox(Surv(time, event) ~ p(fmale, linear=TRUE, tv=T), data=dat), 
    "factor variables")
})

test_that("linear time-varying effect of scalar works with NAs", {
  set.seed(121212)
  eta <- matrix(sin(2 * pi*(1:J)/J) %x% male, nrow=N, ncol=J)
  dat <- simTVSurv(eta, data.frame(z=z, male = male))
  dat$malena <- dat$male; dat$malena[1:5] <- NA
  fit_na <- pcox(Surv(time, event) ~ p(malena, linear=TRUE, tv=T), data=dat)
  expect_gt(cor(coef(fit_na)$value, sin(2*pi*(coef(fit_na)$t/J))), .94)
})  

test_that("multiple linear time-varying effects of scalars work", {
  set.seed(121212)
  eta <- matrix(sin(2*pi*(1:J)/J) %x% z, nrow=N, ncol=J) +
    matrix(cos(2*pi*(1:J)/J) %x% male, nrow=N, ncol=J)
  dat <- simTVSurv(eta, data.frame(z=z, male = male))
  fit <- pcox(Surv(time, event) ~ p(z, linear=TRUE, tv=TRUE) +
      p(male, linear=TRUE, tv=TRUE), data=dat)
  est_z <- coef(fit, 1)
  est_male <- coef(fit, 2)
  expect_gt(cor(est_z$value, sin(2 * pi * est_z$t/J)), .84)
  expect_gt(cor(est_male$value, cos(2 * pi * est_male$t/J)), .94)
  expect_equal(predict(fit), predict(fit, newdata=dat, stimes = dat$time))
  expect_gt(cor(as.vector(predict(fit)), 
    as.vector(eta[, sort(unique(Surv(dat$time, dat$event)[dat$event == 1]))]), 
    use = "pairwise"), .96)
})  

if(FALSE) {
# do these actually work... ? throw rank deficient X warning regardless of "k",
# estimates are quite bad (cor(pred, eta) is .2 for N=500, .45 for N=5000, 
#  ranges are very different), see #43
  # ------------------------------------------------------------------------------
  context("estimate nonlinear time-varying scalar effects")
  
  test_that("nonlinear time-varying effect of scalar works", {
    set.seed(121212)
    eta <- matrix(sin(2 * pi*(1:J)/J) %x% (pnorm(z) - .5), nrow=N, ncol=J)
    dat <- simTVSurv(eta, data.frame(z=z, male = male, fmale=factor(male)))
    fit <- pcox(Surv(time, event) ~ p(z, linear=FALSE, tv=T), data=dat)
    est <- coef(fit)
    pred <- predict(fit)
    expect_gt(cor(as.vector(pred),
      as.vector(eta[, sort(unique(Surv(dat$time, dat$event)[dat$event == 1]))]),
      use = "pairwise"),  .9)
  })
}
