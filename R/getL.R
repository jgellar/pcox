
getL <- function(tind, integration, n.int=NULL) {
  nt <- ncol(tind)
  if (is.null(n.int)) {n.int=rep(nt,nrow(tind))}
  L <- t(sapply(1:nrow(tind), function(i) {
    nt.i <- n.int[i]
    tind.i <- tind[i,1:nt.i]
    L.i <- switch(integration, simpson = {
      ((tind.i[nt.i] - tind.i[1])/nt.i)/3 * c(1, rep(c(4, 
                                                       2), length = nt.i - 2), 1)
    }, trapezoidal = {
      diffs <- diff(tind.i)
      if (length(diffs)>1) {
        0.5 * c(diffs[1], filter(diffs, filter=c(1,1))[-(nt.i-1)],
                diffs[(nt.i-1)])
      } else {
        rep(0.5*diffs,2)
      }
    }, riemann = {
      diffs <- diff(tind.i)
      c(mean(diffs), diffs)
    })
    L.i <- c(L.i, rep(0,nt-nt.i))
  }))
  L
}
