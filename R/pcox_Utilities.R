
#' Deparse a function safely
#' 
#' @param An expression
#' 
#' This utility function was copied from the \code{refund} package.

safeDeparse <- function (expr) {
  ret <- paste(deparse(expr), collapse = "")
  gsub("[[:space:]][[:space:]]+", " ", ret)
}
