library(devtools)
library(survival)
library(mgcv)
library(matlab)
library(fields)

dev_mode()
load_all()

# Function to get the estimate of a historical cox model
getHCEst <- function(sm, idx, coefs, trim=NULL) {
  if (!is.null(trim))
    coefs <- coefs[grepl(trim, names(coefs))]
  J <- length(idx)
  tmat <- matrix(idx, nrow=J, ncol=J)
  smat <- matrix(idx, nrow=J, ncol=J, byrow=TRUE)
  tmat <- as.vector(tmat[lower.tri(tmat, diag=TRUE)])
  smat <- as.vector(smat[lower.tri(smat, diag=TRUE)])
  pmat <- PredictMat(sm, data=data.frame(tmat=tmat, smat=smat, LX=1))
  est  <- matrix(nrow=J, ncol=J, byrow=FALSE)
  est[lower.tri(est,diag=TRUE)] <- as.vector(pmat %*% coefs)
  est
}


# ONE SIMULATION
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
est3 <- getHCEst(sm.out, 1:101, coefs = coef(fit3))

par(mfrow=c(1,2))
image(t(beta1), zlim=c(-6,6), col=jet.colors(64), xaxt="n", yaxt="n", main="True Beta")
axis(1, seq(0,1,length=5), seq(0,100,length=5))
axis(2, seq(0,1,length=5), seq(0,100,length=5))
abline(0,1)
image(t(est3), zlim=c(-6,6), col=jet.colors(64), xaxt="n", yaxt="n", main="Estimate")
axis(1, seq(0,1,length=5), seq(0,100,length=5))
axis(2, seq(0,1,length=5), seq(0,100,length=5))
abline(0,1)


# SOFA Data
data(sofa)
sofa.dead <- sofa[sofa$death,]
fit.sofa1 <- coxph(Surv(los,death) ~ tt(SOFA) + age, data=sofa.dead, na.action=na.pass,
                  tt=tt.func)
sm1 <- sm.out
Jp  <- 50
inc <- 10
est.sofa1 <- getHCEst(sm1, seq(0,Jp,length=101), coef(fit.sofa1), trim="SOFA")

# image.plot(t(est.sofa1), xaxt="n", yaxt="n", zlim=c(-0.1,0.1),
#            main="SOFA:\nMortality", col=jet.colors(64),
#            axis.args=list(at=seq(-.1,.1,by=.05),
#                           labels=c("<-0.10",seq(-.05,.05,by=.05),">0.10")))
image.plot(t(est.sofa1), xaxt="n", yaxt="n", zlim=c(-0.1,0.1),
           main="SOFA:\nMortality", col=jet.colors(64))
image(t(est.sofa1), add=TRUE, zlim=c(-10,-.1), col=jet.colors(64)[1])
image(t(est.sofa1), add=TRUE, zlim=c(.1,1000), col=jet.colors(64)[64])
lbls <- seq(0,Jp,by=inc)
axis(1, lbls/Jp, lbls)
axis(2, lbls/Jp, lbls)
# abline(0,1)

par(mfrow=c(2,5), oma=c(0,0,2,0))
for (i in seq(10,Jp,by=10))
  plot(est.sofa1[i,] ~ c(1:Jp), type="l", main=paste0("t=",i), xlab="")
mtext("Mortality Estimate", side = 3, outer = TRUE)


sofa.dced <- sofa[!sofa$death,]
fit.sofa2 <- coxph(Surv(los,!death) ~ tt(SOFA) + age, data=sofa.dced, na.action=na.pass,
                   tt=tt.func)
sm2 <- sm.out
Jp <- 50
inc <- 10
est.sofa2 <- getHCEst(sm2, seq(0,Jp,length=101), coef(fit.sofa2), trim="SOFA")

image.plot(t(est.sofa2), xaxt="n", yaxt="n", zlim=c(-0.075,.075),
           main="SOFA:\nDischarge", col=jet.colors(64))
image(t(est.sofa2), add=TRUE, zlim=c(-10,-.075), col=jet.colors(64)[1])
image(t(est.sofa2), add=TRUE, zlim=c(.075,1000), col=jet.colors(64)[64])
axis(1, (seq(0,Jp,by=inc)/Jp), seq(0,Jp,by=inc))
axis(2, (seq(0,Jp,by=inc)/Jp), seq(0,Jp,by=inc))
# abline(0,1)

par(mfrow=c(2,5), oma=c(0,0,2,0))
for (i in seq(10,100,by=10))
  plot(est.sofa2[i,] ~ c(1:100), type="l", main=paste0("t=",i), xlab="")
mtext("Hospital Discharge", side = 3, outer = TRUE)


# 2 HEATMAPS
pdf("../Plots/SOFA-HCM-heatmaps50.pdf", width=10, height=5)
par(mfrow=c(1,2), oma=c(0,0,0,2))
image.plot(t(est.sofa1), xaxt="n", yaxt="n", zlim=c(-0.1,0.1),
           main="SOFA:\nMortality", col=jet.colors(64),
           axis.args=list(cex.axis=.75))
image(t(est.sofa1), add=TRUE, zlim=c(-10,-.1), col=jet.colors(64)[1])
image(t(est.sofa1), add=TRUE, zlim=c(.1,1000), col=jet.colors(64)[64])
lbls <- seq(0,Jp,by=inc)
axis(1, lbls/Jp, lbls)
axis(2, lbls/Jp, lbls)
ticks <- (-3:3)*25/1000
image.plot(t(est.sofa2), xaxt="n", yaxt="n", zlim=range(ticks),
           main="SOFA:\nDischarge", col=jet.colors(64),
           axis.args=list(at=ticks, labels=ticks, cex.axis=.75))
image(t(est.sofa2), add=TRUE, zlim=c(-10,-.075), col=jet.colors(64)[1])
image(t(est.sofa2), add=TRUE, zlim=c(.075,1000), col=jet.colors(64)[64])
axis(1, (seq(0,Jp,by=inc)/Jp), seq(0,Jp,by=inc))
axis(2, (seq(0,Jp,by=inc)/Jp), seq(0,Jp,by=inc))
dev.off()


# Bootstrap CIs
B <- 1000
N <- nrow(sofa)
est.sofa1.B = est.sofa2.B <- array(dim = c(101,101,B))
for (b in 1:B) {
  samp <- sample(N, replace=TRUE)
  sofa.b <- sofa[samp,]
  sofa.dead.b <- sofa.b[sofa.b$death,]
  sofa.dced.b <- sofa.b[!sofa.b$death,]
  
  fit1.b <- update(fit.sofa1, data=sofa.dead.b)
  est.sofa1.B[,,b] <- getHCEst(sm.out, seq(0,Jp,length=101), coef(fit1.b), trim="SOFA")
  fit2.b <- update(fit.sofa2, data=sofa.dced.b)
  est.sofa2.B[,,b] <- getHCEst(sm.out, seq(0,Jp,length=101), coef(fit2.b), trim="SOFA")
  
  par(mfrow=c(1,2))
  image.plot(t(est.sofa1.B[,,b]))
  image.plot(t(est.sofa2.B[,,b]))
}

est.sofa1.Bse <- apply(est.sofa1.B, c(1,2), sd)
est.sofa2.Bse <- apply(est.sofa2.B, c(1,2), sd)
est.sofa1.T <- est.sofa1/est.sofa1.Bse
est.sofa2.T <- est.sofa2/est.sofa2.Bse

# Line plots, with SE's
pdf("../Plots/SOFA-HCM-lines50.pdf", width=10, height=8)
par(mfrow=c(2,2))
pal <- jet.colors(11)
plot(est.sofa1[1,1] ~ 0, type="p", pch="o", xlim=c(0,Jp), ylim=c(-.1,.1),
     xlab="s", ylab=expression(beta(s,t)), col=pal[1])
abline(h=0, col="grey")
for (i in 1:10) {
  lines(est.sofa1[10*i+1,] ~ seq(0,Jp,length=101), col=pal[i+1])
}
plot(est.sofa2[1,1] ~ 0, type="p", pch="o", xlim=c(0,Jp), ylim=c(-.1,.1),
     xlab="s", ylab="", col=pal[1])
abline(h=0, col="grey")
for (i in 1:10) {
  lines(est.sofa2[10*i+1,] ~ seq(0,Jp,length=101), col=pal[i+1])
}

plot(est.sofa1.T[1,1] ~ 0, type="p", pch="o", xlim=c(0,Jp), ylim=c(-4,4),
     xlab="s", ylab="T-Statistic", col=pal[1])
abline(h=c(-1.96,0,1.96), col="grey", lty=c("dashed","solid","dashed"))
for (i in 1:10) {
  lines(est.sofa1.T[10*i+1,] ~ seq(0,Jp,length=101), col=pal[i+1])
}
plot(est.sofa2.T[1,1] ~ 0, type="p", pch="o", xlim=c(0,Jp), ylim=c(-4,4),
     xlab="s", ylab="", col=pal[1])
abline(h=c(-1.96,0,1.96), col="grey", lty=c("dashed","solid","dashed"))
for (i in 1:10) {
  lines(est.sofa2.T[10*i+1,] ~ seq(0,Jp,length=101), col=pal[i+1])
}
dev.off()



# Competing risks model
sofa2 <- cbind(id=1:520, sofa)
sofa2 <- rbind(sofa, sofa)
sofa2$eventtype <- c(rep("death", 520), rep("discharge", 520))
sofa2$status <- as.numeric(c(sofa$death, !sofa$death))
fit.cr1 <- coxph(Surv(los,status!=0) ~ strata(eventtype) + tt(SOFA) + age,
                 data=sofa2, na.action=na.pass, tt=tt.func)

