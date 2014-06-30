#' Simulate smooth time-varying (funcitonal) covariates
#' 
#' This function simulates functional covariates, which could also be
#' considered to be a densely-sampled time-varying covariate.
#' The covariates are not truncated,
#' so the function will return a full grid of values that the covariate
#' could take if observed for that long.
#' 
#' @param N sample size (number of rows of \code{X})
#' @param s vector of time points at which the covariate is observed.
#' The length of \code{s}, \code{J}, is used for the number of columns
#' of \code{X}. Defaults to \code{0:100}.
#' 
#' @details This function generates a functional covariate \code{X}
#' from the following model:
#' 
#' @export
#' @author Jonathan Gellar <jgellar1@@jhu.edu>
#' @return An \code{N x J} matrix \code{X}
#' @seealso \code{\link{simTVSurv}}
#' 
genX <- function(N, s=NULL) {
  if (is.null(s)) {
    s  <- (0:100)/100
  }
  u  <- rnorm(N, sd=1)
  v1 <- sapply(1:10, function(k) rnorm(N, sd=2/k))
  v2 <- sapply(1:10, function(k) rnorm(N, sd=2/k))
  X <- sapply(s, function(x) {
    u + rowSums(sapply(1:10, function(k) {
      v1[,k]*sin(2*pi*k*x) + v2[,k]*cos(2*pi*k*x)
    }))
  })
  X
}