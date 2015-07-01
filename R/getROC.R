#' Generate Time-Varying ROC Curve
#' 
#' 
#' @importFrom reshape2 melt
#' @importFrom risksetROC CoxWeights IntegrateAUC
#' @importFrom survival survConcordance.fit
# @export
#' @keywords internal
#' 

getROC <- function(fit, utime=NULL, evaltimes=NULL, plot=FALSE, K=NULL, data=NULL,
                   bhaz=c("none", "breslow"), St=NULL, tmax=NULL) {
  
  # Unique time points
  if (is.null(utime))
    utime <- c(0, rev(unique(fit$y[,1])))
  else if (!(0 %in% utime))
    utime <- c(0, utime)
  
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
    #stop("Not yet implemented!")
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
  tmp <- survConcordance.fit(Surv(Stime, status), lp,
                             as.numeric(factor(Stime, levels=unique(Stime))),
                             NULL)
  if (is.matrix(tmp))
    tmp <- colSums(tmp)
  tmp <- (tmp[1] + tmp[3]/2)/sum(tmp[1:3])
  
  ret <- list(roc=roc, roc.long=roc.long, auc=auc, concordance=tmp)
  #if (!is.null(St)) {
  #  ret$concordance2 <- IntegrateAUC(auc, et1, St, tmax)
  #}
  ret
}


