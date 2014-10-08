getROC <- function(fit, surv, tt, sm, x.var, etimes=NULL) {
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