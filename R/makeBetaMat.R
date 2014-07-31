#' Functions for simulating beta coefficient matrices
#' 
#' @export

makeBetaMat <- function(J=101, bfunc=genBeta1, ...) {
  grid <- do.call("rbind", sapply(0:(J-1), function(t) {
    s <- 0:t
    cbind(s,t)
  }, simplify=FALSE))
  # bvec <- bfunc(grid[,1]/max(grid[,2]), grid[,2]/max(grid[,2]))
  bvec <- bfunc(grid[,1], grid[,2], ...)
  bmat <- matrix(NA,nrow=J,ncol=J)
  for (i in 1:length(bvec)) {
    bmat[grid[i,2]+1,grid[i,1]+1] <- bvec[i]
  }
  bmat
}

genBeta1 <- function(t,Ti,scale=NULL) {
  # Proportion-based, no interaction
  10*ifelse(Ti==0, 0, (t/Ti)-.5)
}

genBeta2 <- function(t,Ti,scale=NULL) {
  # Proportion-based, linear interaction
  if (is.null(scale)) {scale <- max(Ti)}
  10*ifelse(Ti==0, 0, (1-2*Ti/scale)*(0.5-4*(t/Ti-.5)^2))
}

genBeta3 <- function(t,Ti,scale=NULL) {
  # Lag-based, logistic function
  u <- (Ti-t)
  5 - u/10
}

genBeta4 <- function(t,Ti,scale=NULL) {
  # Lag-based, logistic function
  u <- (Ti-t)
  sin(2*pi*Ti/scale)*(5 - u/10)
}

genBeta5 <- function(s,t,lag=3) {
  ifelse((t-s)==lag, 5, 0)
}
