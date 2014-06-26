#' Simulate time-varying survival data
#' 
#' 

simTVSurv <- function (Xdat, etaMat, timeIndex=NULL,
                       eventRandom=NULL, censorRandom=NULL, groupByD=FALSE)
{
  if (!is.data.frame(Xdat) & !is.matrix(Xdat)) {
    stop("Xdat must be a data frame or matrix")
  } else if (!is.data.frame(Xdat)) {
    Xdat <- as.data.frame(Xdat)
  }
  
  numSubjects <- nrow(Xdat)
  if (is.null(timeIndex)) {
    timeIndex <- 1:ncol(etaMat)
  } else if (length(timeIndex) != ncol(etaMat))
    stop("length of timeIndex does not match the number of columns of etaMat")
  maxTime <- timeIndex[length(timeIndex)]
  
  # Generate survival and censoring times
  if (is(eventRandom, "NULL")) 
    eventRandom <- function(n) sample(maxTime, numSubjects, replace = TRUE)
  if (is(eventRandom, "function")) {
    survivalTime <- as.integer(eventRandom(numSubjects))
  } else if (is(eventRandom, "numeric")) {
    if (length(eventRandom) != numSubjects) 
      stop("length of eventRandom is not equal to number of subjects")
    survivalTime <- as.integer(eventRandom)
  } else stop("eventRandom is neither numeric nor function")
  
  if (is(censorRandom, "NULL")) 
    censorRandom <- function(n) sample(maxTime, n, replace = TRUE)
  if (is(censorRandom, "function")) {
    censorTime <- as.integer(censorRandom(numSubjects))
  } else if (is(censorRandom, "numeric")) {
    if (length(censorRandom) != numSubjects) 
      stop("length of censorRandom is not equal to number of subjects")
    censorTime <- as.integer(censorRandom)
  } else stop("censorRandom is neither numeric nor function")
  
  if (min(survivalTime) <= 0) 
    stop("Not all event times are positive")
  notCensored <- ifelse(survivalTime <= apply(cbind(censorTime, 
                                                    maxTime), MARGIN = 1, FUN = min), 1, 0)
  observedTime <- apply(cbind(survivalTime, censorTime, maxTime), 
                        MARGIN = 1, FUN = min)
  
  # Permutation Step
  I <- count <- integer(0)
  if (groupByD) {
    h <- 2 * observedTime + notCensored
    hBins = tabulate(h)
    I <- count <- integer(length(hBins))
    I[h] <- seq(along.with = h)
    I <- I[I > 0]
    count[I] <- hBins[h[I]]
  } else {
    I <- order(observedTime)
    count = rep(1, times = length(I))
  }
  p = .PermuteCovs(observedTime, notCensored, count, 
                                I, etaMat)
  
  # Format output
  if (groupByD) {
    p[order(order(observedTime), sample(numSubjects))]
    J = order(h)
    tuples = cbind(event = notCensored[J], time = observedTime[J], 
                   id.tuples = J)
  } else {
    tuples = cbind(event = notCensored[I], time = observedTime[I],
                   id.tuples = I)
  }
  info = data.frame(cbind(tuples, cov.id = p))
  ordered.info = info[order(info$cov.id), ]
  rownames(ordered.info) <- NULL
  return(.formData(ordered.info, maxTime, numSubjects, Xdat))
}

.PermuteCovs <- function (t, d, count, I, etaMat) 
{
  n = sum(count[I])
  p <- integer(n)
  v <- seq(n)
  ip = 1
  for (k in I) {
    if (d[k]) {
      # Not censored
      J = sample(length(v), size = count[k], replace = FALSE, 
      #           prob = .partialHazards3(t[k], v, etaMat))
      prob = exp(etaMat[v,t[k]]))
    }
    else {
      # Censored
      J = sample(length(v), size = count[k], replace = FALSE, 
                 prob = NULL)
    }
    p[seq(from = ip, along.with = J)] <- v[J]
    ip <- ip + count[k]
    v = v[-J]
  }
  return(p)
}


.formData <- function (ordered.info, m, n, Xdat) 
{
  Xdat <- as.data.frame(sapply(Xdat, function(x) {
    if (is.matrix(x)) {
      return(t(sapply(1:n, function(i) {
        t.i <- ordered.info$time
        c(x[i,1:t.i], rep(NA, m-t.i))
      })))
    } else return(x)
  }, simplify=FALSE))
  data <- cbind(ordered.info[,c("event","time")], Xdat)
  return(data)
}
