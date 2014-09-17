#' Process a time-varying term in a \coe{pcox} formula
#' 
#' 
#' 

tv <- function(term, env=NULL, label=NULL, method="aic", eps=.001,
               basistype="s", bs="tp", ...) {
  # Processing for time-varying functions
  
  
  # Create the appropriate "tt" function to pass to coxph
  tt.func <- function(x,t,...) {
    # Start with code to create the smooth object here
    # Then:
    assign(label, sm, env)
    pterm(sm, method=method, eps=eps)
  }
  
  
}
