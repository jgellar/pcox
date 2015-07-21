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
  any.smooth <- as.logical(length(object$pcox$smooth))
  
  if (any.smooth) {
    names(object$pterms)[!sapply(object$pcox$smoothmap, is.null)] <-
      sapply(object$pcox$smooth, function(sm) sm$label)
  }
  class(object) <- class(object)[-1]
  ret <- summary(object, ...)
  
  if (!is.null(object$pterms)) {
    # Modify conf.int to remove pterms
    idxs <- unlist(object$assign2[!object$pterms])
    ret$conf.int <- ret$conf.int[idxs,,drop=FALSE]
  }
  
  ret$any.smooth <- any.smooth
  ret$class <- class(ret)
  class(ret) <- "summary.pcox"
  ret
}


#' @export
print.summary.pcox <- function(x, ...) {
  class(x) <- x$class
  print(x, ...)
  if ("summary.coxph.penal" %in% class(x) & x$any.smooth) {
    cat("NOTE: smooth pcox coefficients and SEs may be extracted with coef()\n\n")
  }
}


