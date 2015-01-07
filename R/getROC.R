getROC <- function(fit, surv, evaltimes=NULL, plot=FALSE, xval=FALSE,
                   bhaz=c("none", "breslow")) {
  # Get risk set information
  rs <- getRiskSet(surv, evaltimes)
  lp <- fit$linear.predictors
  
  #same results if we limit rs to the predict.time days?  YES!!!!!!!!
  # if (is.null(evaltimes)) -> use info from fit (fit$linear.predictors & fit$y)
  # otherwise, call predict() using evaltimes
  # if evaltimes given:
  #    Calculate linear predictor at each evaltime, for all in risk set, using predict()
  
  p.time <- 34.5
  cw1 <- CoxWeights(marker=as.vector(lp), entry=rs$start, Stime=rs$finish,
             status=rs$newstat, predict.time=p.time)
  ooo <- rs$finish>p.time & rs$start<p.time
  rs.o <- rs[ooo,]
  lp.o <- lp[ooo]
  cw2 <- CoxWeights(marker=as.vector(lp.o), entry=rs.o$start, Stime=rs.o$finish,
                    status=rs.o$newstat, predict.time=p.time)
  
  
  
  
  
  if (is.null(evaltimes)) {
    lp <- fit$linear.predictors
    evaltimes <- unique(rs$finish)
  } else {
    stop("Not yet supported")
    
    environment(tt)$sm.in <- sm
    pt <- tt(x.var[rs$id,], rs$finish)
    environment(tt)$sm.in <- NULL
    pt <- cbind(pt, sofa[rs$id, c("age", "male", "Charlson")])
    pt$male <- as.numeric(pt$male)
    lp <- (as.matrix(pt) - matrix(fit$means, nrow=nrow(pt),
                                  ncol=ncol(pt), byrow=TRUE)
    ) %*% fit$coefficients
  }
  
  # Get ROC curve at each unique time point
  roc <- lapply(evaltimes, function(p.time) {
    CoxWeights(marker=as.vector(lp), entry=rs$start, Stime=rs$finish,
               status=rs$newstat, predict.time=p.time)
  })
  names(roc) <- evaltimes
  roc.long <- do.call("rbind", lapply(1:length(evaltimes), function(i) {
    roc.i <- roc[[i]]
    data.frame(TP=roc.i$TP, FP=roc.i$FP, t=evaltimes[i])
  }))
  auc <- sapply(roc, function(x) x$AUC)
  list(roc=roc, roc.long=roc.long, auc=auc)
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