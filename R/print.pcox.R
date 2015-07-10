#' @export
print.pcox <- function(x, ...) {
  x$call <- x$pcox$call
  names(x$pterms)[(x$pterms==1)] <- sapply(x$pcox$smooth, function(sm) sm$label)
  #names(x$pterms) <- c("A", "B")
  class(x) <- class(x)[-1]
  print(x, ...)
  if (any(x$pterms==1))
    cat("\nNOTE: smooth pcox coefficients may be extracted with coef()\n")
}

#' @export
summary.pcox <- function(object, ...) {
  if (any(object$pterms==1))
    cat("NOTE: smooth pcox coefficients may be extracted with coef()\n\n")
  
  object$call <- object$pcox$call
  names(object$pterms)[(object$pterms==1)] <-
    sapply(object$pcox$smooth, function(sm) sm$label)
  class(object) <- class(object)[-1]
  summary(object, ...)
  
}