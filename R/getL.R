#' Get the weight matrix for a linear functional term
#' 
#' @keywords internal
#' 

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
      } else if (length(diffs)==1) {
        rep(0.5*diffs,2)
      } else {
        1
      }
    }, riemann = {
      diffs <- diff(tind.i)
      c(mean(diffs), diffs)
    })
    L.i <- c(L.i, rep(0,nt-nt.i))
  }))
  L
}

getL2 <- function(tind, integration, n.int=NULL) {
  nt <- ncol(tind)
  if (is.null(n.int)) {n.int=rep(nt,nrow(tind))}
  L <- t(sapply(1:nrow(tind), function(i) {
    nt.i <- n.int[i]
    tind.i <- tind[i,1:nt.i]
    L.i <- switch(integration, simpson = {
      ((tind.i[nt.i] - tind.i[1])/nt.i)/3 * c(1, rep(c(4, 
                                                       2), length = nt.i - 2), 1)
    }, trapezoidal = {
      diffs <- diff(tind.i[1:nt.i])
      if (length(diffs)>1) {
        
        0.5 * c(diffs[1], filter(diffs, filter=c(1,1))[-(nt.i-1)],
                diffs[(nt.i-1)])
      } else if (length(diffs)==1) {
        rep(0.5*diffs,2)
      } else {
        1
      }
    }, riemann = {
      diffs <- diff(tind.i)
      c(mean(diffs), diffs)
    })
    L.i <- c(L.i, rep(0,nt-nt.i))
  }))
  L
}

getL3 <- function(tind, integration, mask=NULL) {
  nt <- ncol(tind)
  if (!is.null(mask)) {
    mask[!mask] <- NA
    tind <- tind*mask
  }
  t(apply(tind, 1, function(tind.i) {
    tvec <- tind.i[!is.na(tind.i)]
    nt <- length(tvec)
    L.i <- switch(integration, simpson = {
      ((tind[nt]-tind[1])/nt)/3 * matrix(c(1, rep(c(4,2), length=nt-2), 1),
                                         nrow=n, ncol=nt, byrow=TRUE)
    }, trapezoidal = {
      diffs <- diff(tvec)
      if (length(diffs)>1) {
        0.5 * c(diffs[1], filter(diffs, filter=c(1,1))[-(nt-1)],
                diffs[(nt-1)])
      } else if (length(diffs)==1) {
        rep(0.5*diffs,2)
      } else {
        1
      }
    }, riemann = {
      diffs <- diff(tvec)
      c(mean(diffs), diffs)
    })
    tind.i[!is.na(tind.i)] <- L.i
    tind.i[is.na(tind.i)] <- 0
    tind.i
  }))
}