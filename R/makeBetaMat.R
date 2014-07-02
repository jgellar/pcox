
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

