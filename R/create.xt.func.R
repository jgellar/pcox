#' Create x-transform function for non-time-varying terms
#' 
#' @keywords internal
#' 

#create.xt.func <- function(limits, linear, basistype, sind, integration,
#                           standardize, domain, basisargs, method, eps) {
create.xt.func <- function(limits, linear, basistype, sind, integration,
                           standardize, s.transform, t.transform,
                           basisargs, method, eps) {
  
  # Pre-processing
  if (!is.null(limits))
    # limits is either "all", "full", or "baseline" - baseline function
    limits <- function(s,t) {s==s} # Will just spit out TRUE's
  smooth.flag <- !(is.null(limits) & linear) # No smooth if its a linear scalar
  smooth  <- NULL
  
  s.transform <- if (is.null(s.transform))
    function(s, t) s
  else if (is.character(s.transform)) {
    
  }
  
  xt.func <- function(x) {
    if (smooth.flag) {
      # Smooth involved: call pcoxTerm
      pcoxTerm(x, limits=limits, linear=linear, tv=FALSE,
               basistype=basistype, sind=sind,
               integration=integration, standardize=standardize,
               s.transform=s.transform, t.transform=t.transform,
               basisargs=basisargs, method=method, eps=eps, smooth=smooth)
    } else {
      # No smooth involved: just return the data (glorified identity function)
      x
    }
  }
  xt.func
}

# If smooth, the return is either (a) a list w/ coxph.penalty and smooth
# objects, or (b) a prediction matrix - if smooth is supplied.
# If not smooth, the return is just the data