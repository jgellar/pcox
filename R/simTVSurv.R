#' Simulate time-varying survival data
#' 
#' This is a modification of the \code{PermAlgo::permalgorithm()}
#' function, which implements the permutation algorithm of Sylvestre
#' and Abrahamowicz to simulate survival data with time-varying
#' coefficients. This version simulates time-varying survival data
#' under more general conditions. In particular, the function can
#' simulate data with time-varying coefficients, time-varying effects,
#' or combinations of the two.
#' 
#' @param etaMat an \code{N x J} matrix contining the linear predictor
#' for each subject at each time index, where \code{N} is the number
#' of subjects and \code{J} is the number of time points.
#' @param Xdat an option \code{N x p} data frame or matrix, containing the
#' covariates used in the model, where \code{p} is the number of
#' covariates. Any covariates supplied will be included in the output
#' dataset. Time-varying covariates are specified as an
#' \code{N x J} matrix using the \code{\link{I}} ("AsIs") function;
#' see Details.
#' @param timeIndex A vector of length \code{J} containing the time
#' indices corresponding to the columns of \code{etaMat}; defaults
#' to \code{1:J}. THIS IS NOT YET TESTED FOR TIME INDICES OTHER THAN
#' \code{1:J} - BEST LEFT AS \code{NULL} FOR NOW.
#' @param eventRandom individual event times. Can either be a vector
#' non-negative integer values, or a random generating function
#' with argument \code{n}. See Details.
#' @param censorRandom individual censoring times, similar to
#' \code{eventRandom}.
#' @param groupByD Logical, should individual assignment of event and
#' censoring times be replaced with group assignments? Defaults to
#' \code{FALSE}. See Details.
#' 
#' @details We assume the following Cox-type regression model:
#' 
#' \eqn{log(h_i(t)) = log(h_0(t)) + \eta_i(t)}
#' 
#' where \eqn{h_i(t)} is the hazard function for subject \eqn{i}
#' at time \eqn{t}, \eqn{h_0(t)} is the baseline hazard function,
#' and \eqn{\eta_i(t)} is the linear predictor matrix, which is
#' a function of the covariates. As opposed to
#' \code{permalgorithm}, which \eqn{\eta_i(t)} by having the
#' user supply a vector of \eqn{\beta} coefficients and assuming
#' a concurrent, non-time-varying effect for time-varying
#' covariates, \code{simTVSurv} accepts \eqn{\eta_i(t)} directly
#' as an argument. Though this approach requires the user to do
#' some more calculations outside of this function, it allows added
#' flexibility by letting the user create the matrix using whatever
#' model they prefer. In particular, this allows the user to
#' specify time-varying effects for both time-fixed and time-varying
#' covariates, and non-concurrent effects of time-varying covariates.
#' 
#' Any covariates included in \code{Xdat} will be added to the output
#' dataset. Any time-varying predictors must first be organized into an
#' \code{N x J} matrix, indicating a value for the covariate for
#' each subject at each \code{timeIndex}. The matrix is then included
#' as a "column" of \code{Xdat} (which must be a data frame), by
#' specifying it as an "AsIs" term with \code{I()} (see examples).
#' When survival/censoring times are generated, time-varying covariates
#' will be "trimmed" into a "time-varying matrix," by replacing
#' any values after the event/censoring time with \code{NA}.
#' 
#' \code{eventRandom} and \code{censorRandom} can each be specified as
#' either a vector of non-negative integer values, or as a random
#' generating function with argument \code{n}. In either case, the
#' values must be smaller or equal to \code{J}. Both
#' default to Uniform[1,\code{J}].
#' \code{groupByD} is an option that, when enabled, increases the
#' compuational efficiency of the algorithm by replaceing the
#' individual assignment of event times and censoring times by
#' grouped assignments. The side effect of this option is that it
#' generates datasets that are, on average, slightly less consistent
#' with the model described by \code{etaMat} than those generated when 
#' \code{groupByD=FALSE}. Still, \code{groupByD=TRUE} may be useful
#' to generate large datasets where \code{J<<N} so that many ties
#' are expected. 
#' 
#' @export
#' @author Jonathan Gellar <jgellar1@@jhu.edu>
#' @return A data frame containing the simulated data, as well as
#' any predictors supplied in \code{Xdat}. Time-varying covariates
#' will be returned as a "time-varying matrix" in "AsIs" format.
#' @seealso \code{\link{permalgorithm}}
#' 

simTVSurv <- function (etaMat, Xdat=NULL, timeIndex=NULL,
                       eventRandom=NULL, censorRandom=NULL, groupByD=FALSE)
{
  N <- nrow(etaMat)
  if (is.null(timeIndex)) {
    timeIndex <- 1:ncol(etaMat)
  } else if (length(timeIndex) != ncol(etaMat))
    stop("length of timeIndex does not match the number of columns of etaMat")
  J <- timeIndex[length(timeIndex)]
  
  if (!is.null(Xdat)) {
    if (!is.data.frame(Xdat) & !is.matrix(Xdat)) {
      stop("Xdat must be a data frame or matrix")
    } else if (!is.data.frame(Xdat)) {
      Xdat <- as.data.frame(Xdat)
    }
    if (nrow(Xdat)!=N)
      stop("etaMat and Xdat must have the same number of rows")
  }
  
  # Generate survival and censoring times
  if (is(eventRandom, "NULL")) 
    eventRandom <- function(n) sample(J, N, replace = TRUE)
  if (is(eventRandom, "function")) {
    survivalTime <- as.integer(eventRandom(N))
  } else if (is(eventRandom, "numeric")) {
    if (length(eventRandom) != N) 
      stop("length of eventRandom is not equal to number of subjects")
    survivalTime <- as.integer(eventRandom)
  } else stop("eventRandom is neither numeric nor function")
  
  if (is(censorRandom, "NULL")) 
    censorRandom <- function(n) sample(J, n, replace = TRUE)
  if (is(censorRandom, "function")) {
    censorTime <- as.integer(censorRandom(N))
  } else if (is(censorRandom, "numeric")) {
    if (length(censorRandom) != N) 
      stop("length of censorRandom is not equal to number of subjects")
    censorTime <- as.integer(censorRandom)
  } else stop("censorRandom is neither numeric nor function")
  
  if (min(survivalTime) <= 0) 
    stop("Not all event times are positive")
  notCensored <- ifelse(survivalTime <= apply(cbind(censorTime, J),
                                              MARGIN = 1, FUN = min), 1, 0)
  observedTime <- apply(cbind(survivalTime, censorTime, J), 
                        MARGIN = 1, FUN = min)
  J <- max(observedTime)
  
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
  p = permuteCovs(observedTime, notCensored, count, I, etaMat)
  
  # Format output
  if (groupByD) {
    p[order(order(observedTime), sample(N))]
    ord = order(h)
    tuples = cbind(event = notCensored[ord], time = observedTime[ord], 
                   id.tuples = ord)
  } else {
    tuples = cbind(event = notCensored[I], time = observedTime[I],
                   id.tuples = I)
  }
  info = data.frame(cbind(tuples, cov.id = p))
  data = info[order(info$cov.id), "event", "time"]
  rownames(ordered.info) <- NULL
  if (!is.null(Xdat))
    data <- addCovData(data, J, N, Xdat)
  return(data)
}

permuteCovs <- function (t, d, count, I, etaMat) 
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


addCovData <- function (data, J, n, Xdat) {
  Xdat <- as.data.frame(sapply(Xdat, function(x) {
    if (is.matrix(x)) {
      if (ncol(x)<J)
        stop("Time-varying covariates must have at least J columns")
      return(I(t(sapply(1:n, function(i) {
        t.i <- data$time
        c(x[i,1:t.i], rep(NA, J-t.i))
      }))))
    } else return(x)
  }, simplify=FALSE))
  data <- cbind(data, Xdat)
  return(data)
}





