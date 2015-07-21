#' @export
print.pcox <- function(x, ...) {
  x$call <- x$pcox$call
  if (length(x$pcox$smoothmap)) {
    names(x$pterms)[!sapply(x$pcox$smoothmap, is.null)] <-
      sapply(x$pcox$smooth, function(sm) sm$label)
  }
  class(x) <- class(x)[-1]
  print(x, ...)
  if (any(x$pterms==1))
    cat("\nNOTE: smooth pcox coefficients may be extracted with coef()\n")
}

#' @export
summary.pcox <- function(object, ...) {
  object$call <- object$pcox$call
  if (length(object$pcox$smoothmap)) {
    names(object$pterms)[!sapply(object$pcox$smoothmap, is.null)] <-
      sapply(object$pcox$smooth, function(sm) sm$label)
  }
  class(object) <- class(object)[-1]
  summ <- summary(object, ...)
  
  if (!is.null(object$pterms)) {
    # Modify conf.int
    idxs <- unlist(object$assign2[!object$pterms])
    summ$conf.int <- summ$conf.int[idxs,,drop=FALSE]
  }
  
  
  # Print summary
  print(summ)
  if (any(object$pterms==1))
    cat("NOTE: smooth pcox coefficients and SEs may be extracted with coef()\n\n")
  invisible(summ)
}
