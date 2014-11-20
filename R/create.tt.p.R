#' Create tt function for a scalar term
#' 
#' @keywords internal
#' 

create.tt.p <- function(limits, linear, tv, basistype, sind,
                        basisargs, method, eps, map,
                        divide.by.t=FALSE, domain=c("s", "s-t", "u"),
                        integration=c("simpson", "trapezoidal", "riemann")) {
  
  smooth <- NULL
  
  # Process limits argument: create appropriate processing function
  conc.fcn <- NULL
  if (is.character(limits)) {
    if (limits %in% c("s=t", "s==t", "t")) {
      #conc.fcn <- function(s,t) s==t
      conc.fcn <- function(s,t) abs(s-t)
      limits <- NULL
    } else if (tolower(limits) %in% c("all", "full")) {
      limits <- function(s,t) {s==s}
    } else if (limits=="s<t") {
      limits <- function(s,t) s<t
    } else if (limits=="s<=t") {
      limits <- function(s,t) s<=t
    } else {
      stop("Unrecognized limits option!")
    }
  } else if (is.numeric(limits)) {
    lag <- limits
    #conc.fcn <- function(s,t) {s==(t-lag)}
    conc.fcn <- function(s,t) abs(s - (t-lag))
    limits <- NULL
    #stop("Numeric limits not currently supported... coming soon")
  }
  
  tt.func <- function(x, t, ...) {
    
    # Turn x matrix back into a data.frame according to map
    if (is.vector(x)) {
      data <- data.frame(x)
    } else {
      data <- as.data.frame(lapply(map, function(map.i) {
        if (length(map.i>1))
          I(x[,map.i])
        else
          x[,map.i]
      }))
    }
    names(data) <- names(map)
    
    # Process concurrent time-varying covariates
    if (!is.null(conc.fcn)) {
      # Reduce data matrices to vectors
      cfidx <- sapply(data, is.matrix)
      data[cfidx] <- lapply(data[cfidx], function(x) {
        if (ncol(x)!=length(sind))
          stop("Mismatch between length of sind and number of columns of data matrix")
        #if (any(sapply(t, function(t.i) min(conc.fcn(sind, t.i)))>1e-2))
        #  stop("Mismatch between event times and sind")
        idxs <- sapply(t, function(t.i) which.min(conc.fcn(sind, t.i)))
        #idxs <- sapply(t, function(t.i) which.min(abs(sind-t.i)))
        sapply(1:nrow(x), function(i) x[i,idxs[i]])
      })
      names(data)[cfidx] <- paste0(names(data)[cfidx], ".t")
    }
    
    if (is.null(limits) & !tv & linear) {
      # No penalized term required: just return data
      as.matrix(data)
    } else {
      # Create coxph.penalty term via pcoxTerm()
      pt <- pcoxTerm(data, limits, linear, tv, basistype, sind, basisargs,
                     method, eps, t)
      env$smooth[[i]] <- pt$smooth
      pt$cpobj
    }
  }
}
