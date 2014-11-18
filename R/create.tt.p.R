#' Create tt function for a scalar term
#' 
#' @keywords internal
#' 

create.tt.p <- function(termname, limits, nonlinear, tv, basistype, sind, basisargs,
                        divide.by.t=FALSE, domain=c("s", "s-t", "u"), dbug=FALSE,
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
    stop("Numeric limits not currently supported... coming soon")
  }
  
  tt.func <- function(data, t, ...) {
    
    if (!is.null(conc.fcn)) {
      # Concurrent time-varying covariate: reduce data matrices to vectors
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
      data
    } else {
      # Create coxph.penalty term via pcoxTerm()
      pt <- pcoxTerm()
      env$smooth[[i]] <- pt$smooth
      pt$cpobj
    }
  }
}


create.tt.p <- function(X, tind=NULL, basistype = c("s", "te", "t2"),
                         tv=TRUE,
                         nonlinear=FALSE,
                         divide.by.t = FALSE,
                         domain = c("s", "s-t", "u"),
                         #domain = c("s", "s-t", "s/t"),
                         dbug=FALSE,
                         integration = c("simpson", "trapezoidal", "riemann"),
                         limits=NULL, splinepars=NULL, ...) {
  
  basistype <- match.arg(basistype)
  domain <- match.arg(domain)
  integration <- match.arg(integration)
  
  if (is.null(limits))
    # Default
    limits <- "s<=t"
  
  lag <- NULL
  if (!is.function(limits)) {
    if (is.numeric(limits)) {
      # limits specifies a lag for inclusion
      lag <- limits
      limits <- function(s, t) {s<=t & s>=(t-lag)}
    } else if (limits == "s<t") {
      limits <- function(s, t) {s < t}
    } else if (limits == "s<=t") {
      limits <- function(s, t) {(s <= t)}
    } else {
      stop("supplied <limits> argument unknown")
    }
  }
  sm <- NULL
  
  tt.func <- function(x.var,t.var,...) {
    n <- nrow(x.var)
    J <- ncol(x.var)
    tmat <- matrix(t.var, nrow=n, ncol=J)
    smat <- matrix(1:J, nrow=n, ncol=J, byrow=TRUE)
    mask <- t(outer(smat[1,], tmat[,1], limits))
    
    if (domain=="s-t") {
      # Align according to s-t
      smat <- smat-tmat
    } else if (domain=="u") {
      stop("Not currently supported")
      
      # Domain-standardized
      # Problem when t==min(smat)?
      # Need different L matrix? NO (I think....)
      minS <- min(smat)
      smat[t.var!=minS,] <- (smat[t.var!=minS,]-minS)/(tmat[t.var!=minS,]-minS)
      smat[t.var==minS,] <- 0.5
      
      #smat <- (smat-min(smat))/(tmat-min(smat))
    }
    
    L <- getL3(smat, integration=integration, mask=mask)
    x.var[is.na(x.var)] <- 0
    rownames(x.var) <- NULL
    LX <- L*x.var
    if (divide.by.t)
      LX <- if (is.null(lag)) {
        LX/matrix(t.var, nrow=n, ncol=J)
      } else {
        LX <- LX/matrix(sapply(t.var, min, lag), nrow=n, ncol=J)
      }
    
    call <- list(as.symbol(basistype), quote(smat))
    evaldat <- data.frame()
    evaldat$smat <- smat
    
    # If time-varying, add tmat
    if (tv) {
      call <- c(call, quote(tmat))
      evaldat$tmat <- tmat
    }
    
    # If nonlinear, add xmat within s(), otherwise include in "by" variable
    if (nonlinear) {
      call <- c(call, x.var, by=quote(L))
      evaldat$x.var <- x.var
      evaldat$L <- L
    } else {
      LX <- L*x.var
      call <- c(call, by=quote(LX))
      evaldat$LX <- LX
    }
    
    # Add formals....
    
    
    
    if (is.null(sm)) {
      # Create new smooth object
      sm <- smoothCon(eval(as.call(call)), data=evaldat,
                      knots=NULL, absorb.cons=TRUE)[[1]]
      
      # Assign it to the env environment
      assign(smooth[[index]], sm, env)
      
      # Create coxph.penalty term
      pterm(sm, method=method, eps=eps)
    } else {
      # Call PredictMat() using existing smooth object on new data
      PredictMat(sm, data = evaldat)
    }
    
    #pdat <- data.frame(tmat=I(tmat), smat=I(smat), LX=I(LX))
    
    #if (is.null(sm)) {
    #  sm <- if (is.null(k)) {
    #    smoothCon(s(tmat, smat, by=LX), data=pdat,
    #              knots=NULL, absorb.cons=TRUE)[[1]]
    #  } else {
    #    smoothCon(s(tmat, smat, by=LX, k=k), data=pdat,
    #              knots=NULL, absorb.cons=TRUE)[[1]]
    #  }
    #  # NEED TO GET RID OF <<-
    #  sm.out <<- sm
    #  pterm(sm, method=method, eps=eps)
    #} else {
    #  PredictMat(sm, data = pdat)
    #}
  }
  
  if (dbug) {
    debug(tt.func)
  } else if (isdebugged(tt.func))
    undebug(tt.func)
  tt.func  
}

