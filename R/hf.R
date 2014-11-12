
hf <- function(X, ...) {
  tt.func <- create.tt.hf(X, ...)
  attr(X, "tt") <- tt.func
  X
}

