#' Create tt function for a scalar term
#' 
#' @keywords internal
#' 

create.tt.sm <- function(additive = TRUE, tv = FALSE, basistype = c("s", "te", "t2"),
                         dbug=FALSE, ...) {
  
  basistype <- match.arg(basistype)
  if (!additive & !tv)
    stop("Nothing to smooth over for one of the s terms!")
  dots <- list(...)
  sm <- NULL
  
  tt.func <- function(x.var, t.var=NULL, ...) {
    evaldat <- data.frame(x.var=x.var)
    if (!is.null(t.var)) {
      evaldat$t.var <- t.var
    } else if (tv==TRUE) {
      stop("times must be supplied for time-varying terms")
    }
    
    if (is.null(sm)) {
      # Set up appropriate call to construct basis
      #s2 <- mgcv::s
      #newcall <- list(ifelse(basistype=="s", quote(s2), as.symbol(basistype)))
      newcall <- list(as.symbol(basistype))
      
      # If additive, add x to the smooth
      if (additive) {
        newcall <- c(newcall, quote(x.var))
        pfx <- 1
      } else {
        pfx <- x.var
      }
      
      # If tv, add time to smooth
      if (tv) {
        newcall <- c(newcall, quote(t.var))
      }
      
      # Add dots...????
      newcall <- c(newcall, dots)
      
      # Create new smooth object
      sm <- smoothCon(eval(as.call(newcall)), data=evaldat,
                      knots=NULL, absorb.cons=TRUE)[[1]]
      
      # Assign it to the env environment
      env$smooth[[index]] <- sm
      #assign(smooth[[index]], sm, env)
      
      # Create coxph.penalty term (pre-multiply by x.var if not additive)
      pfx * pterm(sm, method=method, eps=eps)
    } else {
      # Call PredictMat() using existing smooth object on new data
      PredictMat(sm, data = evaldat)
    }    
  }
  
  if (dbug) {
    debug(tt.func)
  }
  tt.func
}


create.tt.sm2 <- function(..., additive = TRUE, tv = FALSE,
                         basistype = c("s", "te", "t2"),
                         dbug=FALSE, k=NULL, fx=NULL, bs=NULL, m=NULL,
                         by=NULL, xt=NULL, id=NULL, sp=NULL) {
  
  basistype <- match.arg(basistype)
  if (!additive & !tv)
    stop("Nothing to smooth over for one of the s terms!")
  sm <- NULL
  
  tt.func <- function(x.var, t.var=NULL, ...) {
    evaldat <- data.frame(x.var=x.var)
    if (!is.null(t.var)) {
      evaldat$t.var <- t.var
    } else if (tv==TRUE) {
      stop("times must be supplied for time-varying terms")
    }
    
    if (is.null(sm)) {
      # Set up appropriate call to construct basis
      s2 <- mgcv::s
      newcall <- list(ifelse(basistype=="s", quote(s2), as.symbol(basistype)))
      
      if (additive) {
        newcall <- c(newcall, quote(x.var))
        pfx <- 1
      } else {
        pfx <- x.var
      }
      
      if (tv) {
        newcall <- c(newcall, quote(t.var))
      }
      
      # Add dots...????
      #newcall <- c(newcall, list(...))
      
      # Create new smooth object
      sm <- smoothCon(eval(as.call(newcall)), data=evaldat,
                      knots=NULL, absorb.cons=TRUE)[[1]]
      
      # Assign it to the env environment
      env$smooth[[index]] <- sm
      #assign(smooth[[index]], sm, env)
      
      # Create coxph.penalty term (pre-multiply by x.var if not additive)
      pfx * pterm(sm, method=method, eps=eps)
    } else {
      # Call PredictMat() using existing smooth object on new data
      PredictMat(sm, data = evaldat)
    }    
  }
  
  if (dbug) {
    debug(tt.func)
  }
  
  tt.func
}

