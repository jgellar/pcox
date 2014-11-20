# pcox::sm

sm <- function(x, additive = TRUE, tv = FALSE, basistype = c("s", "te", "t2"),
               dbug=FALSE, ...) {
  tt.func <- if (additive | tv) {
    create.tt.sm(additive = additive, tv = tv, basistype = basistype, dbug=dbug, ...)
  } else {
    NULL
  }
  list(x=x, tt=NULL)
}


sm2 <- function(..., additive = TRUE, tv = FALSE, basistype = c("s", "te", "t2"),
               dbug=FALSE) {
  frmls <- formals(match.fun(basistype))
  dots  <- list(...)
  infrmls <- names(dots) %in% names(frmls)
  x    <- dots[!infrmls]
  dots <- dots[ infrmls]
  
  x2 <- as.data.frame(lapply(x, function(y) ifelse(is.null(dim(y)), y, I(y))))
  
  x2 <- as.data.frame(ifelse(sapply(x, function(y) is.null(dim(y))),
                             x, I(x)))
                     
  
  tt.func <- if (additive | tv) {
    arglist <- c(additive=additive, tv=tv, basistype=basistype, dbug=dbug, dots)
    do.call(create.tt.sm, arglist)
  } else NULL
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

