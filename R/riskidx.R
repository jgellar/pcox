riskidx <- function (stimes, status, etimes=NULL) {
  if (length(stimes)!=length(status))
    stop("stimes and status must have same length")
  
  if (min(stimes[status])<=0) {
    stimes[status & stimes<=0] <- 1e-04
    warning("Survival times <= 0 replaced with 1e-04")
  }
  
  utimes <- if (is.null(etimes)) {
    unique(stimes[status!=0])
  } else {
    unique(c(stimes[(stimes <= max(etimes)) & status!=0], etimes))
  }
  utimes <- utimes[order(utimes)]
  nt = length(utimes)
  nc = length(stimes)
  t.evaluate = c(0, utimes)
  
  as.data.frame(do.call("rbind", sapply(1:nt, function(i) {
    start   <- t.evaluate[i]
    finish  <- t.evaluate[i+1]
    keep    <- stimes>=finish
    newstat <- ifelse(stimes[keep]==finish, status[keep], 0)
    id <- which(keep)
    cbind(start, finish, newstat, id)
  })))
}