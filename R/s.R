# pcox::s

s <- function(x, ...) {
  tt.func <- create.tt.s(x, ...)
  list(x=x, tt=tt.func)
  #attr(x, "tt") <- tt.func
  #x
}
