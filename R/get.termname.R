#' Accessory function to extract the term name from a tt or xt function
#' @keywords internal
#' 

get.termname <- function(func, varnames=NULL) {
  env <- environment(func)
  
  if (is.null(varnames))
    # get varnames from tt function map
    varnames <- names(env$map)
  
  # Process concurrent TVC's
  if (!is.null(env$conc.fcn)) {
    # Add .t to appropriate variable names
    cfidx <- sapply(map, length) > 1
    varnames[cfidx] <- paste0(varnames[idx], ".t")
  }
  
  # xt functions won't have a $tv variable
  if (is.null(env$tv))
    env$tv <- FALSE
  
  # Build term name
  if (is.null(env$limits) & !(env$tv) & env$linear) {
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
      if (env$linear)
        byvar <- paste0(varnames[1], ".LX")
      else {
        inside <- c(inside, varnames)
        byvar <- paste0(varnames[1], ".L")
      }
    }    
    
    paste0(env$basistype, "(", paste(inside, collapse = ", "),
                 ifelse(is.null(byvar), "", paste0(", by = ", byvar)),
                 ")")
  }
}