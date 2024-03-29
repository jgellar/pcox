# coef.pcox: old version

oldcoef.pcox <- function(object, raw=FALSE, term=NULL, inds=NULL, se=FALSE) {
    if (raw) {
      # Raw coefficients only
      if (se) {
        return(list(coefficients=object$coefficients,
                    se1=sqrt(diag(object$var)),
                    se2=sqrt(diag(object$var2))))
      } else {
        return(object$coefficients)
      }
    } else {
      term.labels <- attr(object$pcox$terms, "term.labels")
      if (is.null(term)) {
        # All terms (as a list)
        res <- lapply(term.labels, function(x) {
          ttype <- a
        })
        
      } else {
        # One term
        if (is.character(term)) {
          term <- pmatch(term, term.labels)
        }
        if (!is.numeric(term)) stop("Error: Unrecognized term")
        
        
        #res <- acmCoef(object, ??, inds, se)
        
        
        
        
        smap.i <- object$pcox$smoothmap[[i]]
        if (is.null(smap.i)) {
          # No smooth involved - just return coefficient(s)
          data.frame()
        } else {
          # Smooth term: call appropriate mgcv:::plot function
          smooth.i <- object$smooth[[]]
        }
        
      }
    }  
}


# 
# coef.acm <- function(object, raw=FALSE, term=NULL, inds=NULL, se=FALSE) {
#   if (raw) {
#     # Raw coefficients only
#     if (se) {
#       return(list(coefficients=object$coefficients,
#                   se1=sqrt(diag(object$var)),
#                   se2=sqrt(diag(object$var2))))
#     } else {
#       return(object$coefficients)
#     }
#   } else {
#     term.labels <- attr(object$acm$terms, "term.labels")
#     if (is.null(term)) {
#       # All terms (as a list)
#       res <- lapply(term.labels, function(x) {
#         ttype <- a
#       })
#       
#     } else {
#       # One term
#       if (is.character(term)) {
#         term <- pmatch(term, term.labels)
#       }
#       if (!is.numeric(term)) stop("Error: Unrecognized term")
#       
#       
#       #res <- acmCoef(object, ??, inds, se)
#       
#     }
#   }
# }
# 
# acmCoef <- function(coef, sm, inds=NULL, var1=NULL, var2=NULL) {
#   
# }
# 




