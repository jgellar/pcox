#' Accessory function that allows s(), te(), and t2() in coxph term names
#' @keywords internal
#' 

# Returns the variable that matches the call, instead of calling the
# corresponding function. Useful for allowing special term names.
f_override <- function(...) {
  callstring <- deparse(match.call(), width.cutoff = 500L)
  if (length(callstring) > 1) callstring <- paste0(callstring)
  get(callstring, parent.frame())
}