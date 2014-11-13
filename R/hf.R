
hf <- function(X, ...) {
  tt.func <- create.tt.hf(X, ...)
  list(x=X, tt=tt.func)
  #attr(X, "tt") <- tt.func
  #X
}

