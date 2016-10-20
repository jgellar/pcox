library(pcox)
source(system.file("tests/setup-testthat.R", package = "pcox"))

test_that("parametric model works", {
  # h(x) = .6*x + .75*male 
  eta <- matrix(.6*x + .75*male, nrow=N, ncol=J)
  dat <- simTVSurv(eta, data.frame(x=x, male=male))
  fit <- pcox(Surv(time, event) ~ x + male, data=dat)
  fit_coxph <- survival::coxph(Surv(time, event) ~ x + male, data=dat)
  expect_equivalent(coef(fit), coef(fit_coxph))
  expect_equivalent(predict(fit),  predict(fit_coxph))
  expect_equivalent(predict(fit),  predict(fit, newdata=dat))
  
  dat_NA <- dat
  dat_MNA
  
})

test_that("smooth term works", {
  # h(x) = sin(x) + .75*male 
  eta <- matrix(sin(x) + .75*male, nrow=N, ncol=J)
  dat <- simTVSurv(eta, data.frame(x=x, male=male))
  fit <- pcox(Surv(time, event) ~ p(x, linear=FALSE, k=8) + male, data=dat)
  
  
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

