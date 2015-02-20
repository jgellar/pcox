#' Create tt function for a time-varying terms
#' 
#' @keywords internal
#' 


create.tt.func <- function(limits, linear, tv, basistype, sind, integration,
                           standardize, s.transform, t.transform,
                           basisargs, method, eps, map) {
  
  # Initialize: no smooth object created yet
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
    if (limits<0)
      stop("limits must be non-negative if using it to specify lagged 
           time-varying covariates")
    lag <- limits
    conc.fcn <- function(s,t) abs(s - (t-lag))
    limits <- NULL # Treat the lagged covariate as a scalar
  } else if (!is.function(limits)) {
    stop("Unrecognized input for limits argument")
  }
  
  # Process optional t transformation
  t.transform <- if (is.null(t.transform)) {
    if (!is.null(s.transform) & basistype=="s") {
      if (s.transform=="s/t")
        # Defaults to "tmaxmin"
        function(t, tmax, tmin) (t-tmin)/(tmax-tmin)
      else NULL            # Defaults to no transform
    } else NULL            # Defaults to no transform
  } else if (is.character(t.transform)) {
    if (t.transform=="t") NULL
    else if (t.transform=="tmax") function(t, tmax) t/tmax
    else if (t.transform=="tmaxmin") function(t, tmax, tmin)
      (t-tmin)/(tmax-tmin)
    else stop("Unrecognized t transformation")
  } else if (!is.function(t.transform))
    stop("Unrecognized t tranformation: must be a function or a
         recognized transformation string")
  else if (length(formals(t.transform))>3 | length(formals(t.transform))<1)
    stop("t.transform must have 1, 2, or 3 arguments")
  
  # Process optional s transformation
  s.transform <- if (is.null(s.transform))
    # Defaults to no transform
    NULL
  else if (is.character(s.transform)) {
    if (s.transform=="s") NULL
    else if (s.transform=="s-t") function(s,t) s-t
    else if (s.transform=="s/t") function(s,t) {
      smin <- min(s)
      ifelse(t==smin, 0.5, (s-smin)/(t-smin))
    }
    else stop("Unrecognized s transformation")
  } else if (!is.function(s.transform))
    stop("Unrecognized s tranformation: must be a function or a
         recognized transformation string")
  else if (length(formals(s.transform))>2 | length(formals(s.transform))<1)
    stop("s.transform must have 1 or 2 arguments")
  
  tt.func <- function(x, t, ...) {
    
    # Turn x matrix back into a data.frame according to map
    if (is.vector(x)) {
      data <- data.frame(x)
    } else {
      rownames(x) <- NULL
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
          stop("Mismatch between length of sind and number of columns 
               of data matrix")
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
      #   (i.e., the most basic concurrent TVC - no extras)
      as.matrix(data)
    } else {
      # Create pcoxTerm
      pt <- pcoxTerm(data, limits, linear, tv, basistype, sind, integration,
                     standardize, s.transform, t.transform,
                     basisargs, method, eps, smooth, t)
      if (is.list(pt)) {
        # tt.func is being called within pcox/coxph, to create a new term:
        #   Assign the smooth, and return cpobj
        env$smooth[[index]] <- pt$smooth
        pt$cpobj
      } else {
        # tt.func is called for a pre-existing term:
        #   Return the prediction matrix
        pt
      }
    }
  }
  # Return
  tt.func
}
