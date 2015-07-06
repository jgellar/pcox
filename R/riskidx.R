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

#' @useDynLib pcox coxcount1
getRiskSet <- function (Y, evaltimes=NULL) {
  if (is.null(evaltimes)) {
    sorted <- order(-Y[,1], Y[,2])
    newstrat <- rep.int(0L, nrow(Y))
    newstrat[1] <- 1L
    if (storage.mode(Y) != "double") 
      storage.mode(Y) <- "double"
    counts <- .Call(coxcount1, Y[sorted,], as.integer(newstrat))
    rs <- data.frame(start=rep(c(counts$time[-1],0), counts$nrisk),
                     finish=rep(counts$time, counts$nrisk),
                     newstat=counts$status, id=sorted[counts$index])
    #rs2 <- arrange(rs2, finish, id)
  } else {
    stop("Not yet supported!")
  }  
  rs
}