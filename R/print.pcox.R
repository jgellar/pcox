#' @export
print.pcox <- function(object, ...) {
  object$call <- object$pcox$call
  #names(object$pterms) <- c("A", "B")
  class(object) <- class(object)[-1]
  print(object, ...)
}

#' @export
summary.pcox <- function(object, ...) {
  object$call <- object$pcox$call
  class(object) <- class(object)[-1]
  summary(object, ...)
}