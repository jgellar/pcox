# pcox::s

s <- function(x, ...) {
  tt.func <- create.tt.s(x, ...)
  list(x=x, tt=tt.func)
}

s3 <- function(..., additive = TRUE, tv = FALSE,
               basistype = c("s", "te", "t2"), dbug=FALSE,                
               k=-1, fx=FALSE, bs="tp", m=NA, by=NA, xt=NULL, id=NULL, sp=NULL) {
  x <- data.frame(...)
  tt.func <- create.tt.s(x, additive = additive, tv = tv,
                         basistype = basistype, dbug=dbug,                
                         k=k, fx=fx, bs=bs, m=m, by=by, xt=st, id=id, sp=sp)
  list(x=x, tt=tt.func)
}

