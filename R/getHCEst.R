#' Extract a historical Cox estimate from a fitted model
#' 
#' @export
#' 

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