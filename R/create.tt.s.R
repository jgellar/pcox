
# tt function for smooth scalar term

create.tt.s <- function(x, additive = TRUE, tv = FALSE,
                        basistype = c("s", "te", "t2"),
                        dbug=FALSE, ...) {
  
  basistype <- match.arg(basistype)
  if (!additive & !tv)
    stop("Nothing to smooth over for one of the s terms!")
  sm <- NULL
  
  tt.func <- function(x.var, t.var, ...) {
    call <- list(ifelse(basistype=="s", mgcv::s, match.fun(basistype)))
    evaldat <- data.frame(x.var=x.var, t.var=t.var)
    
    pfx <- 1
    if (additive) {
      call <- c(call, quote(x.var))
    } else {
      pfx <- x.var
    }
    
    if (tv) {
      call <- c(call, quote(t.var))
    }
    
    # Add dots...????
    #call <- c(call, list(...))
    
    if (is.null(sm)) {
      # Create new smooth object
      sm <- smoothCon(eval(as.call(call)), data=evaldat,
                      knots=NULL, absorb.cons=TRUE)[[1]]
      
      # Assign it to the env environment
      assign(smooth[[index]], sm, env)
      
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



