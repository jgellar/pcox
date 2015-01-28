#' Accessory function to extract the term name from a tt function
#' @keywords internal
#' 

get.ttname <- function(tt) {
  env <- environment(tt)
  
  varnames <- names(env$map)
  if (!is.null(env$conc.fcn)) {
    # Add .t to appropriate variable names
    cfidx <- sapply(map, length) > 1
    varnames[cfidx] <- paste0(varnames[idx], ".t")
  }
  
  if (is.null(env$limits) & !env$tv & env$linear) {
    # Basic concurrent TVC - no basistype call
    varnames[1]
  } else {
    # A smooth object will be created: build up term label
    inside <- vector("character", length=0)
    if (is.null(env$limits)) {
      # Smooth of scalar variable(s)
      if (env$tv)
        inside <- "t"
      if (!env$linear) {
        inside <- c(inside, varnames)
        byvar <- NULL
      }
      else
        byvar <- varnames
    } else {
      # Functional Term
      inside <- paste0(varnames[1], ".smat")
      if (env$tv)
        inside <- c(inside, paste0(varnames[1], ".tmat"))
      if (!env$linear) {
        inside <- c(inside, varnames)
        byvar <- paste0(varnames[1], ".L")
      } else
        byvar <- paste0(varnames[1], ".LX")
    }    
    
    nm <- paste0(env$basistype, "(", paste(inside, collapse = ", "), ")")
    ifelse(is.null(byvar), nm, paste0(nm, ":", byvar))
  }
}