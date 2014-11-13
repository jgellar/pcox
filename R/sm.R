# pcox::sm

sm <- function(x, additive = TRUE, tv = FALSE, basistype = c("s", "te", "t2"),
               dbug=FALSE, ...) {
  tt.func <- create.tt.sm(x, additive = additive, tv = tv, basistype = basistype,
                          dbug=dbug, ...)
  list(x=x, tt=tt.func)
}

# sm3 <- function(..., additive = TRUE, tv = FALSE,
#                basistype = c("s", "te", "t2"), dbug=FALSE,                
#                k=-1, fx=FALSE, bs="tp", m=NA, by=NA, xt=NULL, id=NULL, sp=NULL) {
#   x <- data.frame(...)
#   tt.func <- create.tt.s(x, additive = additive, tv = tv,
#                          basistype = basistype, dbug=dbug,                
#                          k=k, fx=fx, bs=bs, m=m, by=by, xt=st, id=id, sp=sp)
#   list(x=x, tt=tt.func)
# }
# 

