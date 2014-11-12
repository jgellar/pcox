#' Extract a historical Cox estimate from a fitted model
#' 
#' @export
#' 

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
  pmat <- PredictMat(sm, data=data.frame(tmat=tmat, smat=smat, LX=1))
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
