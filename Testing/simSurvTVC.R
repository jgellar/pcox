
# Generate random survival data based on a PH model with
# TVCs or TVEs


Y <- ceiling(rweibull(N, .75, 600/gamma(1+1/.75)))
C <- rep(365*2, N)
delta <- Y<C
Y[!delta] <- C[!delta]


# \eta_i(t) = 5 * x_i(t)


# PermAlgo
n=500
m=365
Xmat=matrix(ncol=3, nrow=n*m)
TDhist <- function(m){
  start <- round(runif(1,1,m),0) # individual start date
  duration <-  7 + 7*rpois(1,3) # in weeks
  dose <-  round(runif(1,0,10),1) 
  vec <- c(rep(0, start-1), rep(dose, duration))
  while (length(vec)<=m){
    intermission <- 21 + 7*rpois(1,3) # in weeks
    duration <-  7 + 7*rpois(1,3) # in weeks
    dose <-  round(runif(1,0,10),1)
    vec <- append(vec, c(rep(0, intermission), rep(dose, duration)))}
  return(vec[1:m])}
Xmat[,1] <- rep(rbinom(n, 1, 0.3), each=m)
Xmat[,2] <- do.call("c", lapply(1:n, function(i) TDhist(m)))
Xmat[,3] <- do.call("c", lapply(1:n, function(i) TDhist(m)))

eventRandom <- round(rexp(n, 0.012)+1,0)
censorRandom <- round(runif(n, 1,870),0)
data <- permalgorithm(n, m, Xmat, XmatNames=c("sex", "Drug1", "Drug2"),
                      eventRandom = eventRandom, censorRandom=censorRandom,
                      betas=c(log(2), log(1.04), log(0.99)), groupByD=FALSE )

fit.pa <- coxph(Surv(Start,Stop,Event) ~ sex + Drug1 + Drug2, data)



# Now with my simulated X's: TVC
X <- genX(N, s=seq(0,10,length=365*2))
J <- ncol(X)
Xmat <- matrix(nrow=N*J, ncol=2)
Xmat[,1] <- rep(rbinom(N,1,.5), each=J)
Xmat[,2] <- as.vector(t(X))
eventRandom <- round(rweibull(N, .75, 600/gamma(1+1/.75)))+1
censorRandom <- round(runif(N, 1, 1000))
data <- permalgorithm(N, J, Xmat, XmatNames=c("sex", "SOFA"),
                      eventRandom=eventRandom, censorRandom=censorRandom,
                      betas=c(log(2), log(1.5)))
fit <- coxph(Surv(Start,Stop,Event) ~ sex + SOFA, data)


# permalgorithm2: TVE
t.ind <- 1:J
betas <- rbind(log(2), 3*sin(2*pi*t.ind/J))
Xmat[,2] <- rep(rnorm(N), each=J)
data <- permalgorithm2(N, J, Xmat, XmatNames=c("sex", "SOFA"),
                       eventRandom=eventRandom, censorRandom=censorRandom,
                       betas=betas)
data2 <- do.call("rbind",by(data, data$Id, function(x) {
  x[nrow(x),]
}, simplify=TRUE))
fit <- timecox(Surv(Stop,Event) ~ const(sex) + SOFA, data2)
tmp.t <- fit$cum[,1]
tmp.x <- fit$cum[,3]
est <- diff(predict(gam(tmp.x~s(tmp.t))))/diff(tmp.t)
plot(est ~ tmp.t[-length(tmp.t)], type="l")
lines(betas[2,], col="red")

par(mfrow=c(1,2))
plot(fit)
lines(betas[2,], col="red")
tmpFit <- gam(tmp.x ~ s(tmp.t))


# Try again - large N
N = 5000
X <- genX(N, s=seq(0,10,length=365*2))
J <- ncol(X)
Xmat <- matrix(nrow=N*J, ncol=2)
Xmat[,1] <- rep(rbinom(N,1,.5), each=J)
Xmat[,2] <- as.vector(t(X))
eventRandom <- round(rweibull(N, .75, 600/gamma(1+1/.75)))+1
censorRandom <- round(runif(N, 1, 1000))
data <- permalgorithm2(N, J, Xmat, XmatNames=c("sex", "SOFA"),
                      eventRandom=eventRandom, censorRandom=censorRandom,
                      betas=betas)
data2 <- do.call("rbind",by(data, data$Id, function(x) {
  x[nrow(x),]
}, simplify=TRUE))
fit <- timecox(Surv(Stop,Event) ~ const(sex) + SOFA, data2)
tmp.t <- fit$cum[,1]
tmp.x <- fit$cum[,3]
est <- diff(predict(gam(tmp.x~s(tmp.t))))/diff(tmp.t)
plot(est ~ tmp.t[-length(tmp.t)], type="l", ylim=c(-6,3))
lines(betas[2,], col="red")
data.big <- data2
fit.big <- fit

# Now try with tt()
N <- 500
t.ind <- 1:J
betas <- rbind(log(2), 3*sin(2*pi*t.ind/J))
Xmat <- matrix(nrow=N*J, ncol=2)
Xmat[,1] <- rep(rbinom(N,1,.5), each=J)
Xmat[,2] <- rep(rnorm(N), each=J)
eventRandom <- round(rweibull(N, .75, 600/gamma(1+1/.75)))+1
censorRandom <- round(runif(N, 1, 1000))
data.tve <- permalgorithm2(N, J, Xmat, XmatNames=c("sex", "SOFA"),
                       eventRandom=eventRandom, censorRandom=censorRandom,
                       betas=betas)
data.tve <- do.call("rbind",by(data.tve, data.tve$Id, function(x) {
  x[nrow(x),]
}, simplify=TRUE))
fit1 <- timecox(Surv(Stop,Event) ~ const(sex) + SOFA, data.tve)
tmp.t <- fit1$cum[,1]
tmp.x <- fit1$cum[,3]
est1 <- diff(predict(gam(tmp.x~s(tmp.t))))/diff(tmp.t)
plot(est1 ~ tmp.t[-length(tmp.t)], type="l")
lines(betas[2,], col="red")

fit2 <- coxph(Surv(Stop,Event) ~ sex + tt(SOFA), data=data,
              tt=function(x,t,...) x*pspline(t))
est2 <- predict(fit2, type="terms",
                newdata=data.frame())
fit3 <- coxph(Surv(Stop,Event) ~ sex + tt(SOFA), data=data,
              tt=function(x,t,...) x*s.cox(t))
# tvec from inside s.cox....
sm <- smoothCon(s(tvec), data=data.frame(tvec=tvec),
                knots=NULL, absorb.cons=TRUE)[[1]]
pmat <- PredictMat(sm, data=data.frame(tvec=1:730))
est3 <- as.vector(pmat %*% coef(fit3)[-1])
lines(as.vector(est3) ~ c(1:730), col="blue")

fit4 <- aareg(Surv(Stop,Event) ~ sex + SOFA, data=data)
plot(fit4$times, cumsum(fit4$coef[,2]), type="l")
est4 <- cumsum(fit4$coef[,2])
  


out.a <- aalen(  Surv(time,status==9)~const(age)+const(sex)+
                   const(diabetes)+chf+vf,
                 data=sTRACE,max.time=7,n.sim=100)
out.c <- timecox(Surv(time,status==9)~const(age)+const(sex)+
                   const(diabetes)+chf+vf,
                 data=sTRACE,max.time=7,n.sim=100)
out.s <- coxph(Surv(time,status==9) ~ age+sex+diabetes+
                 tt(chf) + tt(vf), data=sTRACE,
               tt=function(x,t,...) x*s.cox(t))
sm <- smoothCon(s(tvec), data=data.frame(tvec=tvec.cc),
                knots=NULL, absorb.cons=TRUE)[[1]]
pmat <- PredictMat(sm, data=data.frame(tvec=seq(min(sTRACE$time),
                                                    max(sTRACE$time),
                                                    length=200)))
est.s1 <- pmat %*% coef(out.s)[4:12]
est.s2 <- pmat %*% coef(out.s)[13:21]

tmp.t <- out.a$cum[,1]
est.a1 <- diff(predict(gam(c(out.a$cum[,3])~s(tmp.t))))/diff(tmp.t)
est.a2 <- diff(predict(gam(c(out.a$cum[,4])~s(tmp.t))))/diff(tmp.t)

tmp.t <- out.c$cum[,1]
est.c1 <- diff(predict(gam(c(out.c$cum[,3])~s(tmp.t))))/diff(tmp.t)
est.c2 <- diff(predict(gam(c(out.c$cum[,4])~s(tmp.t))))/diff(tmp.t)

plot(est.a1 ~ out.a$cum[-1,1], type="l", ylim=range(est.c1))
lines(est.c1 ~ out.c$cum[-1,1], col=2)
lines(est.s1 ~ seq(min(sTRACE$time), max(sTRACE$time), length=200), col=3)

plot(est.a2 ~ out.a$cum[-1,1], type="l", ylim=c(-2,3))
lines(est.c2 ~ out.c$cum[-1,1], col=2)
lines(est.s2 ~ seq(min(sTRACE$time), max(sTRACE$time), length=200), col=3)


# TVCs
N <- 200
X <- genX(N)
J <- ncol(X)
Xmat <- matrix(nrow=N*J, ncol=2)
Xmat[,1] <- rep(rbinom(N,1,.5), each=J)
Xmat[,2] <- as.vector(t(X))
eventRandom <- round(rweibull(N, .75, 100/gamma(1+1/.75)))+1
censorRandom <- round(runif(N, 1, 200))
data.tvc <- permalgorithm(N, J, Xmat, XmatNames=c("sex", "SOFA"),
                       eventRandom=eventRandom, censorRandom=censorRandom,
                       betas=c(log(2), log(3)))
SOFA <- do.call("rbind", tapply(data.tvc$SOFA, data.tvc$Id, function(x) {
  c(x, rep(NA, J-length(x)))
}))
data.tvc <- do.call("rbind",by(data.tvc, data.tvc$Id, function(x) {
  x[nrow(x),]
}, simplify=TRUE))
data.tvc$SOFA <- I(SOFA)

fit.tvc <- coxph(Surv(Stop,Event) ~ sex + tt(SOFA), data=data.tvc, 
                 tt=tt.func, na.action=na.pass)
fit.tve <- coxph(Surv(Stop,Event) ~ sex + tt(SOFA), data=data.tve,
              tt=function(x,t,...) x*s.cox(t))

sm$first.para <- 2
sm$last.para <- 31
fit.tvc$smooth <- list(sm)

tmat.pre <- matrix(1:J, nrow=J, ncol=J)
smat.pre <- matrix(1:J, nrow=J, ncol=J, byrow=TRUE)
LX.pre  <- matrix(nrow=J,ncol=J)
LX.pre[lower.tri(LX.pre,diag=TRUE)] <- 1
dat.pre <- data.frame(smat=I(smat.pre), tmat=I(tmat.pre), LX=I(LX.pre))
bhat <- getEst(mod=fit.tvc, data=dat.pre)




# TVC model based on simulated data


