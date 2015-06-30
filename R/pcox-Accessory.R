#' \code{getCall} method for \code{pcox} class
#' @keywords internal
getCall.pcox <- function(x) x$pcox$call

#' Handle Missing Values for \code{pcox} Fits
#' 
#' This function differs from \code{na.omit} in the way it handles data frames
#' that contain matrix elements. \code{na.omit} will remove an observation
#' if any of the elements in that row of the matrix are missing, whereas this
#' function only removes the observation if the entire row is missing.
#' @keywords internal
na.omit_pcox <- function(object, ...) {
  n <- length(object)
  omit <- logical(nrow(object))
  vars <- seq_len(n)
  for (j in vars) {
    x <- object[[j]]
    if (!is.atomic(x)) 
      next
    x <- is.na(x)
    d <- dim(x)
    if (is.null(d) || length(d) != 2L) 
      omit <- omit | x
    else
      omit <- omit | apply(x, 1, all)
  }
  xx <- object[!omit, , drop = FALSE]
  if (any(omit > 0L)) {
    temp <- setNames(seq(omit)[omit], attr(object, "row.names")[omit])
    attr(temp, "class") <- "omit"
    attr(xx, "na.action") <- temp
  }
  xx
}

#' Accessory function to extract the term name from a tt or xt function
#' @keywords internal
get.termname <- function(func, varnames=NULL) {
  env <- environment(func)
  
  if (is.null(varnames))
    # get varnames from tt function map
    varnames <- names(env$map)
  
  # Process concurrent TVC's
  if (!is.null(env$conc.fcn)) {
    # Add .t to appropriate variable names
    cfidx <- sapply(env$map, function(x) length(x) > 1)
    varnames[cfidx] <- paste0(varnames[cfidx], ".t")
  }
  
  # xt functions won't have a $tv variable
  if (is.null(env$tv))
    env$tv <- FALSE
  
  # Create term name
  if (is.null(env$limits) & !(env$tv) & env$linear) {
    # Basic concurrent TVC - no basistype call, just return varname[1]
    varnames[1]
  } else {
    # A smooth object will be created: build up term name
    inside <- vector("character", length=0)
    if (is.null(env$limits)) {
      # Smooth of scalar variable(s)
      if (env$tv)
        inside <- "t"
      if (!env$linear) {
        inside <- c(inside, varnames)
        byvar <- NULL
      }
      else
        byvar <- varnames
    } else {
      # Functional Term
      inside <- paste0(varnames[1], ".smat")
      if (env$tv)
        inside <- c(inside, paste0(varnames[1], ".tmat"))
      if (env$linear)
        byvar <- paste0(varnames[1], ".LX")
      else {
        inside <- c(inside, varnames)
        byvar <- paste0(varnames[1], ".L")
      }
    }
    
    # Paste term components and return
    paste0(env$basistype, "(", paste(inside, collapse = ", "),
           ifelse(is.null(byvar), "", paste0(", by = ", byvar)), ")")
  }
}

#' Accessory function that allows s(), te(), and t2() in coxph term names
#' 
#' Returns the variable that matches the call, instead of calling the
#' corresponding function. Useful for allowing special term names.
#' @keywords internal
f_override <- function(...) {
  callstring <- deparse(match.call(), width.cutoff = 500L)
  if (length(callstring) > 1) callstring <- paste0(callstring)
  get(callstring, parent.frame())
}

#' Deparse a function safely
#' 
#' This utility function was copied from the \code{refund} package.
#' @keywords internal
safeDeparse <- function (expr) {
  ret <- paste(deparse(expr), collapse = "")
  gsub("[[:space:]][[:space:]]+", " ", ret)
}

#' Convert a list to a data frame
#' 
#' This function was copied from the \code{refund} package.
#' @keywords internal
list2df <- function (l) {
  nrows <- sapply(l, function(x) nrow(as.matrix(x)))
  stopifnot(length(unique(nrows)) == 1)
  ret <- data.frame(rep(NA, nrows[1]))
  for (i in 1:length(l)) ret[[i]] <- l[[i]]
  names(ret) <- names(l)
  return(ret)
}

#' Expands a call
#' 
#' This function was copied from the \code{refund} package.
#' @keywords internal
expand.call <- function (definition = NULL, call = sys.call(sys.parent(1)), 
          expand.dots = TRUE) {
  call <- match.call(definition, call, expand.dots)
  ans <- as.list(call)
  frmls <- formals(safeDeparse(ans[[1]]))
  frmls <- frmls[!sapply(frmls, is.symbol)]
  add <- which(!(names(frmls) %in% names(ans)))
  return(as.call(c(ans, frmls[add])))
}