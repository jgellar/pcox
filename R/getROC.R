getROC <- function(fit, evaltimes=NULL, plot=FALSE, K=NULL, data=NULL,
                   bhaz=c("none", "breslow")) {
  
  # Unique time points
  utime <- c(0, rev(unique(fit$y[,1])))
  
  # Get lp, Stime, entry, and status
  if (is.null(K)) {
    # No cross-validation
    lp <- fit$linear.predictors
    Stime <- fit$y[,1]
    entry  <- sapply(Stime, function(x)
      utime[which(utime==x)-1])
    status <- as.logical(fit$y[,2])
  } else {
    # K-fold cross-validation
    if (is.null(data))
      stop("data required for cross-validation")
    stop("Not yet implemented!")
    n <- nrow(data)
    folds <- sample(rep(1:K, ceiling(n/K))[1:n])
    pre.mat <- matrix(nrow=n, ncol=(length(utime)-1))
    colnames(pre.mat) <- utime[-1]
    cat("Cross validating model")
    for (k in 1:K) {
      cat(".")
      dat.k <- data[folds!=k,]
      fit.k <- update(fit, data=dat.k)
      dat.k.pre <- data[folds==k,]
      pre.k <- predict(fit.k, newdata=dat.k.pre, ptimes=utime[-1],
                        stimes=data$los[folds==k])
      pre.mat[folds==k, which(colnames(pre.mat) %in% colnames(pre.k))] <- pre.k
    }
    print(" Cross-validation complete")
    pre.long <- melt(pre.mat, varnames = c("Subject", "Stime"), 
                     value.name = "lp", na.rm=TRUE)
    lp <- pre.long$lp
    Stime <- pre.long$Stime
    entry  <- sapply(Stime, function(x)
      utime[which(utime==x)-1])
    finalstatus <- fit$pcox$surv[pre.long$Subject, 2]
    finaltime <- fit$pcox$surv[pre.long$Subject, 1]
    status <- finalstatus & (finaltime==Stime)
  }
  
  if (is.null(evaltimes)) evaltimes <- utime
  
  roc <- lapply(evaltimes, function(t) {
    status.t <- Stime==t & status
    CoxWeights(marker = lp, Stime = Stime, status = status.t,
               predict.time = t, entry = entry)
  })
  names(roc) <- evaltimes
  
  roc.long <- do.call("rbind", lapply(1:length(evaltimes), function(i) {
    roc.i <- roc[[i]]
    data.frame(TP=roc.i$TP, FP=roc.i$FP, t=evaltimes[i])
  }))
  auc <- sapply(roc, function(x) x$AUC)
  list(roc=roc, roc.long=roc.long, auc=auc)
  
  
#   # Get risk set information
#   rs <- getRiskSet(surv, evaltimes)
#   lp <- fit$linear.predictors
#   
#   #same results if we limit rs to the predict.time days?  YES!!!!!!!!
#   # if (is.null(evaltimes)) -> use info from fit (fit$linear.predictors & fit$y)
#   # otherwise, call predict() using evaltimes
#   # if evaltimes given:
#   #    Calculate linear predictor at each evaltime, for all in risk set, using predict()
#   
#   p.time <- 34.5
#   cw1 <- CoxWeights(marker=as.vector(lp), entry=rs$start, Stime=rs$finish,
#              status=rs$newstat, predict.time=p.time)
#   ooo <- rs$finish>p.time & rs$start<p.time
#   rs.o <- rs[ooo,]
#   lp.o <- lp[ooo]
#   cw2 <- CoxWeights(marker=as.vector(lp.o), entry=rs.o$start, Stime=rs.o$finish,
#                     status=rs.o$newstat, predict.time=p.time)
#   
#   if (is.null(evaltimes)) {
#     lp <- fit$linear.predictors
#     evaltimes <- unique(rs$finish)
#   } else {
#     stop("Not yet supported")
#     
#     environment(tt)$sm.in <- sm
#     pt <- tt(x.var[rs$id,], rs$finish)
#     environment(tt)$sm.in <- NULL
#     pt <- cbind(pt, sofa[rs$id, c("age", "male", "Charlson")])
#     pt$male <- as.numeric(pt$male)
#     lp <- (as.matrix(pt) - matrix(fit$means, nrow=nrow(pt),
#                                   ncol=ncol(pt), byrow=TRUE)
#     ) %*% fit$coefficients
#   }
#   
#   # Get ROC curve at each unique time point
#   roc <- lapply(evaltimes, function(p.time) {
#     CoxWeights(marker=as.vector(lp), entry=rs$start, Stime=rs$finish,
#                status=rs$newstat, predict.time=p.time)
#   })
#   names(roc) <- evaltimes
#   roc.long <- do.call("rbind", lapply(1:length(evaltimes), function(i) {
#     roc.i <- roc[[i]]
#     data.frame(TP=roc.i$TP, FP=roc.i$FP, t=evaltimes[i])
#   }))
#   auc <- sapply(roc, function(x) x$AUC)
#   list(roc=roc, roc.long=roc.long, auc=auc)
}

getROC.old <- function(fit, Y, tt, sm, x.var, etimes=NULL) {
  rs <- riskidx(surv[,1], surv[,2], etimes = etimes)
  environment(tt)$sm.in <- sm
  pt <- tt(x.var[rs$id,], rs$finish)
  environment(tt)$sm.in <- NULL
  #pt <- create.tt.func(divide.by.t=TRUE, sm.in=sm)(x.var[rs$id,],
  #                                                 rs$finish)
  pt <- cbind(pt, sofa[rs$id, c("age", "male", "Charlson")])
  pt$male <- as.numeric(pt$male)
  lp <- (as.matrix(pt) - matrix(fit$means, nrow=nrow(pt),
                                ncol=ncol(pt), byrow=TRUE)
  ) %*% fit$coefficients
  t.eval <- unique(rs$finish)
  roc <- lapply(t.eval, function(t) {
    risksetROC2(rs$finish, entry=rs$start, status=rs$newstat,
                marker=as.vector(lp), predict.time=t, plot=F)
  })
  names(roc) <- t.eval
  roc.long <- do.call("rbind", lapply(1:length(t.eval), function(i) {
    roc.i <- roc[[i]]
    data.frame(TP=roc.i$TP, FP=roc.i$FP, t=t.eval[i])
  }))
  auc <- sapply(roc, function(x) x$AUC)
  list(roc=roc, roc.long=roc.long, auc=auc)
}


getROC0 <- function(fit, surv, x.var, etimes=NULL) {
  rs <- riskidx(surv[,1], surv[,2], etimes = etimes)
  t <- rs$finish
  pt <- apply(rs, 1, function(x) x.var[x[4], x[2]])
  pt <- cbind(pt, sofa[rs$id, c("age", "male", "Charlson")])
  pt$male <- as.numeric(pt$male)
  lp <- (as.matrix(pt) - matrix(fit$means, nrow=nrow(pt),
                                ncol=ncol(pt), byrow=TRUE)
  ) %*% fit$coefficients
  t.eval <- unique(rs$finish)
  roc <- lapply(t.eval, function(t) {
    risksetROC2(rs$finish, entry=rs$start, status=rs$newstat,
                marker=as.vector(lp), predict.time=t, plot=F)
  })
  names(roc) <- t.eval
  roc.long <- do.call("rbind", lapply(1:length(t.eval), function(i) {
    roc.i <- roc[[i]]
    data.frame(TP=roc.i$TP, FP=roc.i$FP, t=t.eval[i])
  }))
  auc <- sapply(roc, function(x) x$AUC)
  list(roc=roc, roc.long=roc.long, auc=auc)
}