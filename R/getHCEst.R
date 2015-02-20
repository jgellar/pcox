#' Extract a historical Cox estimate from a fitted model
#' 
#' @export
#' 

getHCEst2 <- function(fit, term=NULL, idx=NULL, mask=NULL) {
  if (is.null(term)) {
    term <- which(fit$pcox$termtype=="hf")[1]
  }
  
  if (is.null(idx))
    idx <- 1:max(fit$y[,1])
  else if (length(idx)==1)
    idx <- 1:idx
  
  sm <- fit$pcox$smooth[[term]][[1]]
  coefs <- fit$coefficients[sm$first.para:sm$last.para]
  env <- environment(fit$pcox$t.funcs[[term]])
  
  J <- length(idx)
  tmat <- matrix(idx, nrow=J, ncol=J)
  smat <- matrix(idx, nrow=J, ncol=J, byrow=TRUE)
  if (is.null(mask))
    mask <- t(outer(smat[1,], tmat[,1], env$limits))
  smat <- as.vector(smat[mask])
  tmat <- as.vector(tmat[mask])
  
  s.transform <- env$s.transform
  t.transform <- env$t.transform
  if (!is.null(s.transform))
    smat <- s.transform(smat, tmat)  
  if (!is.null(t.transform))
    tmat <- t.transform(tmat, max(fit$y[,1]), min(fit$y[,1]))
  
  df <- data.frame(smat=smat, tmat=tmat, LX=1)
  names(df) <- c(sm$term, sm$by)
  pmat <- PredictMat(sm, data=df)
  
  est  <- matrix(nrow=J, ncol=J)
  est[mask] <- as.vector(pmat %*% coefs)
  est
}


# getHCEst2 <- function(fit, term=NULL, idx=NULL, mask=NULL) {
#   if (is.null(term)) {
#     term <- which(fit$pcox$termtype=="hf")[1]
#   }
#   
#   if (is.null(idx))
#     idx <- 1:max(fit$y[,1])
#   else if (length(idx)==1)
#     idx <- 1:idx
#   
#   sm <- fit$pcox$smooth[[term]][[1]]
#   coefs <- fit$coefficients[sm$first.para:sm$last.para]
#   
#   J <- length(idx)
#   if (is.null(mask)) mask <- lower.tri(matrix(nrow=J, ncol=J), diag=TRUE)
#   
#   tmat <- as.vector(matrix(idx, nrow=J, ncol=J)[mask])
#   smat <- as.vector(matrix(idx, nrow=J, ncol=J, byrow=TRUE)[mask])
#   
#   df <- data.frame(smat=smat, tmat=tmat, LX=1)
#   names(df) <- c(sm$term, sm$by)
#   pmat <- PredictMat(sm, data=df)
#   
#   est  <- matrix(nrow=J, ncol=J)
#   est[mask] <- as.vector(pmat %*% coefs)
#   est
# }


getHCEst <- function(sm, idx, coefs, trim=NULL, mask=NULL) {
  if (!is.null(trim))
    coefs <- coefs[grepl(trim, names(coefs))]
  J <- length(idx)
  tmat <- matrix(idx, nrow=J, ncol=J)
  smat <- matrix(idx, nrow=J, ncol=J, byrow=TRUE)
  
  tmat <- if (is.null(mask)) {
    # Default is the lower triangle
    as.vector(tmat[lower.tri(tmat, diag=TRUE)])
  } else {
    # If mask is supplied, use the mask
    as.vector(tmat[mask])
  }
  
  smat <- if (is.null(mask)) {
    as.vector(smat[lower.tri(smat, diag=TRUE)])
  } else {
    as.vector(smat[mask])
  }
  df <- data.frame(smat=smat, tmat=tmat, LX=1)
  names(df) <- c(sm$term, sm$by)
  pmat <- PredictMat(sm, data=df)
  est  <- matrix(nrow=J, ncol=J, byrow=FALSE)
  
  if (is.null(mask)) {
    est[lower.tri(est,diag=TRUE)] <- as.vector(pmat %*% coefs)
  } else {
    est[mask] <- as.vector(pmat %*% coefs)
  }
  est
}


getHCEst.smt <- function(sm, idx, lag, coefs, trim=NULL, mask=NULL, rescale=TRUE) {
  if (!is.null(trim))
    coefs <- coefs[grepl(trim, names(coefs))]
  
  # Set up matrices
  J <- length(idx)
  smt.idx <- seq(-lag,0,by=.5)
  K <- length(smt.idx)
  tmat <- matrix(idx, nrow=J, ncol=K)
  smtmat <- matrix(smt.idx, nrow=J, ncol=K, byrow=TRUE)
  
  # Mask out the lower right corner
  if (is.null(mask))
    mask <- -1*smtmat <= tmat
  tmat <- as.vector(tmat[mask])
  smtmat <- as.vector(smtmat[mask])
  if (rescale) {
    smtmat <- (smtmat-min(smtmat))/(max(smtmat)-min(smtmat))
    tmat <- (tmat-min(tmat))/(max(tmat)-min(tmat))
  }
  
  # Make predictions and organize as matrix
  pmat <- PredictMat(sm, data=data.frame(tmat=tmat, smtmat=smtmat, LX=1))
  est  <- matrix(nrow=J, ncol=K, byrow=FALSE)
  est[mask] <- as.vector(pmat %*% coefs)
  est
}
