# Simulation 1 for Historical Cox Model


N.vec <- seq(100, 200, 500)
m.vec <- c("AIC", "cAIC", "EPIC")

for (r in 1:R) {
  for (i.N in 1:length(N.vec)) {
    N <- N.vec[i.N]
    X <- genX(N)
    
    for (b in 1:length(beta.list)) {
      beta <- beta.list[[b]]
      eta  <- sapply(1:J, function(j) {
        X[,1:j,drop=F] %*% beta1[j,1:j]
        # X[,1:j,drop=F] %*% beta1[j,1:j] / j
      })
      data <- simTVSurv(eta, Xdat=data.frame(X=I(X)))
      
      for (m in 1:length(method.vec)) {
        method <- m.vec[m]
        
        fit <- coxph()
      }
    }
  }
}



J <- 101
beta1 <- makeBetaMat(J, genBeta5, lag=5)

# SCENARIO 3: Historical Functional Terms for TVC's
N <- 500
J <- 101
Xdat <- data.frame(X1=I(genX(N, s=seq(0,1,length=J))), X2=rnorm(N))
beta1 <- makeBetaMat(J, genBeta1)
eta <- sapply(1:J, function(j) {
  1.3*Xdat[[2]] + Xdat[[1]][,1:j,drop=F] %*% beta1[j,1:j] / j
})
data3 <- simTVSurv(eta, Xdat=Xdat)
fit3 <- coxph(Surv(time,event) ~ tt(X1) + X2, data=data3, na.action=na.pass,
              tt=tt.func)
