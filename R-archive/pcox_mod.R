#' Fit a Penalized Cox Regression Model
#' 
#' This function fits a penalized Cox regression model using penalized splines 
#' to model the smooth coefficients. It is a wrapper for \code{survival::coxph}
#' and \code{coxme::coxme} to fit Cox models with smooth covariate effects for
#' (possibly censored) time-to-event data, of the general form \eqn{log[h_i(t)]
#' = log[h_0(t)] + f_1(z_i, t) + f_2(X_i(s), t)} where \eqn{h_i(t)} is the
#' hazard function for subject i, \eqn{h_0(t)} is the baseline hazard function,
#' \eqn{z_i} are scalar covariates, and \eqn{X_i(s)} are functional covariates.
#' 
#' @param formula a formula with special terms as for \code{\link[mgcv]{gam}},
#'   with additional special terms allowed. See Details.
#' @param method method used to optimize the smoothing parameter(s). See
#'   Details.
#' @param eps convergence level for the criterion indicated by \code{method}.
#' @param knots optional list containing the user-specified knot values for
#'   basis construction, as in \code{\link[mgcv]{gam}}.
#' @param x should the design matrix be included in the return object?
#'   Defaults to \code{TRUE}.
#' @param ...   additional arguments for \code{\link[survival]{coxph}} or
#'   \code{\link[coxme]{coxme}}.
#' @details This routine is a wrapper for either \code{\link[survival]{coxph}}
#'   or \code{\link[coxme]{coxme}}, which fit penalized survival models by
#'   maximizing the penalized partial likelihood. The fitting routine is
#'   dependent on the \code{method} chosen for optimizing the smoothing
#'   parameter: ML or REML optimization is performing by \code{coxme}, whereas
#'   likelihood-based criteria (AIC, AIC_c, or EPIC) are implemented by
#'   \code{coxph}. \code{coxph} is also used to implement the penalized survival
#'   model with fixed value of the smoothing parameter (\code{method}="fixed",
#'   and a \code{theta} arguemnt is required) or degrees of freedom
#'   (\code{method}="df", and a \code{df} argument is required).
#'   
#'   Smooth (penalized) terms are indicated in the model formula as follows: 
#'   \enumerate{ \item Nonlinear (and possibly multivariate) effects of scalar
#'   covariates that do not vary with the time index \eqn{t} (i.e., f(z_{i1}))
#'   are specified using either \code{\link{s}}(), or \code{\link{te}}(), as in
#'   \code{\link[mgcv]{gam}}. \item Nonlinear (and possibly multivariate)
#'   effects of scalar covariates that vary with the time index \eqn{t} (i.e.,
#'   f(z_{i1},t)) are specified using ???. \item Linear and nonlinear effects of
#'   non-concurrent functional predictors that do not vary with the time index
#'   \eqn{t} (i.e., \eqn{\int X_i(s)\beta(s)ds} or \eqn{\int\beta(X_i(s), s)
#'   ds}) are specified using \code{\link[refund]{lf}}, 
#'   \code{\link[refund]{af}}, or \code{lf.vd} from the package 
#'   \code{refund}. \item Linear and nonlinear effects of non-concurrent
#'   functional predictors that vary with the time index \eqn{t} are specified
#'   using \code{\link[refund]{ff}} or \code{\link[refund]{sff}} from the
#'   package \code{refund}, respectively. Note that the outcome "function" for
#'   these "function-on-function" terms is the hazard function. \item Concurrent
#'   effects of functional covariates are not yet allowed, but will be soon! }
#'   
#'   Note that we refer to a functional predictor as "concurrent" if it is
#'   measured over the domain \eqn{t}, the same domain as the outcome event is
#'   measured across. Another term for this type of predictor is a
#'   densely-measured time-varying covariate. These types of covariates require
#'   special attention....
#'   
# @import mgcv refund coxme survival
#' @importFrom mgcv gam gam.fit s te t2
#' @importFrom survival coxph Surv
#' @export
#' @author Jonathan Gellar <jgellar1@@jhu.edu> and Fabian Scheipl
#' @return a fitted \code{pcox} object. This is either a \code{coxph} or 
#'   \code{coxme} object with   additional information in the \code{pcox} entry.
#' @references Gellar, Jonathan E., Colantuoni, Elizabeth, Needham, Dale M., and
#'   Crainiceanu, Ciprian M. (2014). "Cox Regression Models with Functional 
#'   Covariates for Survival Data." Johns Hopkins University, Dept. of 
#'   Biostatistics Working Papers. Working Paper 264. 
#'   \url{http://biostats.bepress.com/jhubiostat/paper264}.
#' @seealso \code{mgcv}'s \code{\link[mgcv]{smooth.terms}} for details of 
#'   \code{mgcv} syntax and available spline bases and penalties; the related 
#'   \code{\link[refund]{pffr}} and \code{\link[refund]{fgam}} from 
#'   \code{refund}.
#'   
pcox <- function(formula, data,
                 method=c("aic","caic","epic"),
                 #method=c("aic","caic","epic","reml", "ml", "fixed", "df"),
                 eps=1e-6, knots=NULL, x=TRUE, ...) {
  # Preliminaries...
  call <- match.call()
  method <- match.arg(method)
  dots <- list(...)
  if (!("iter.max" %in% names(dots)))  dots$iter.max  <- 100
  if (!("outer.max" %in% names(dots))) dots$outer.max <- 50
  
  fitter <- ifelse(method %in% c("reml","ml"), "coxme", "coxph")
  #if (length(dots)) {
  #  validDots <- ifelse(fitter=="coxme", names(formals(coxme)),
  #                      names(formals(coxph)))
  #  notUsed <- names(dots)[!(names(dots) %in% validDots)]
  #  if (length(notUsed)) 
  #    warning("Arguments <", paste(notUsed, collapse = ", "), 
  #            "> supplied but not used.")
  #}  
  
  # Organize terms
  tf <- terms.formula(formula, specials = c("p", "bf", "hf", "cf",
                                            "strata", "cluster", "frailty"))
  trmstrings <- attr(tf, "term.labels")
  terms <- sapply(trmstrings, function(trm) as.call(parse(text = trm))[[1]], 
                  simplify = FALSE)
  frmlenv  <- environment(formula)
  specials <- attr(tf, "specials")
  where.p  <- specials$p  - 1
  where.bf <- specials$bf - 1
  where.hf <- specials$hf - 1
  where.cf <- specials$cf - 1
  where.pen <- c(where.p, where.bf, where.hf, where.cf)
  where.par <- if (length(trmstrings)) {
    which(!(1:length(trmstrings) %in% where.pen))
  } else {
    numeric(0)
  }
  
  # Organize response and environment
  responsename <- attr(tf, "variables")[2][[1]]
  newfrml <- paste(safeDeparse(responsename), "~", sep = "")
  newfrmlenv <- new.env()
  evalenv <- if ("data" %in% names(call)) 
    eval.parent(call$data)
  else NULL
  #else parent.frame() # IS THIS OK?
  
  surv <- eval(responsename, envir = evalenv)
  #surv <- eval(responsename, envir = evalenv, enclos = frmlenv)
  nobs <- nrow(surv)
  fitter <- as.symbol(fitter)
  if (is.call(responsename)) {
    # responsename is a call to Surv, assign its arguments to newfrmlenv
    sapply(as.list(responsename[-1]), function(x) {
      assign(x = deparse(x),
             value = eval(x, envir=evalenv, enclos=frmlenv),
             envir = newfrmlenv)
      invisible(NULL)
    })
  } else {
    # assign responsename (the survival object) directly to newfrmlenv
    assign(x = deparse(responsename),
           value = eval(responsename, envir = evalenv, enclos = frmlenv),
           envir = newfrmlenv)
  }
  
  # Set up variables for new formula and smooth objects
  newtrmstrings <- attr(tf, "term.labels")
  t.funcs = varmap = paramap = smooth = smoothdata  <-
    vector(mode = "list", length = length(trmstrings))
  t.types <- rep(NA, length=length(trmstrings))
  #smooth <- vector("list", length=length(newtrmstrings))
  
  # Override s(), te(), and t2() in newfrmlenv so they can be used in
  # term names. This forces s(x), etc. to return the variable "s(x)" instead
  # of calling mgcv::s(). Has to be in the parent of the formula environment
  # so that list2df(newfrmlenv), etc. below still work:
  assign("s",  f_override, envir=parent.env(newfrmlenv))
  assign("te", f_override, envir=parent.env(newfrmlenv))
  assign("t2", f_override, envir=parent.env(newfrmlenv))
  
  
  #################
  # Process Terms #
  #################
  
  # Penalized terms
  if (length(where.pen)) {
    for (i in where.pen) {
      # Evaluate term (which is a call to p() - directly or inderectly)
      terms[[i]]$method <- method
      terms[[i]]$eps <- eps
      trm <- eval(terms[[i]], envir=evalenv)
      #trm <- eval(terms[[i]], envir=evalenv, enclos=frmlenv)
      
      # Extract, modify, and save transformation function
      is.tt <- if (!is.null(trm$tt)) TRUE else if (!is.null(trm$xt)) FALSE else
        stop("Error: shouldn't get here - something's wrong!")
      tf.i <- if (is.tt) trm$tt else trm$xt
      environment(tf.i)$env   <- environment()
      environment(tf.i)$index <- i
      t.funcs[[i]] <- tf.i
      t.types[i] <- ifelse(is.tt, "tt", "xt")
      
      # Process separately for xt/tt
      if (!is.tt) {
        # xt function: execute now
        trm.i <- tf.i(trm$x)
        varmap[[i]] <- names(trm$x)
        nm <- if (class(trm.i)=="coxph.penalty")
          get.termname(tf.i, names(trm$x)) else names(trm.i)
        assign(x=nm, trm.i, envir=newfrmlenv)
        newtrmstrings[i] <- nm
      } else {
        # tt function: execute later (within coxph)
        varmap[[i]] <- names(environment(tf.i)$map)
        nm <- get.termname(tf.i)
        assign(x=nm, trm$x, envir=newfrmlenv)
        newtrmstrings[i] <- paste0("tt(",nm,")")
      }
    }
  }
  
  # Parametric terms
  if (length(where.par)) {
    for (i in where.par) {
      term.i <- terms[i]
      nms <- if (is.call(term.i[[1]]))
        # term is a call: extract variable names
        sapply(as.list(term.i[[1]])[-1], function(cl) as.character(cl) )
      else
        # term is a name, just use it
        as.character(term.i[[1]])
        
      varmap[[i]] <- nms
      sapply(nms, function(nm) {
        v <- eval(as.name(nm), envir = evalenv, enclos = frmlenv)
        stopifnot(length(v) == nobs)
        assign(x = nm, value = v, envir = newfrmlenv)
        invisible(NULL)
      })
    }
  }
  
  
  # Insert code into coxph to remove NAs from model frame AFTER tt terms are processed.
  # This is necessary especially for lagged concurrent effects of TVC's,
  # and possibly for historical effects (if atypical limits are used)
  suppressMessages(
    trace(coxph,
          at=which(sapply(as.list(body(coxph)), function(x)
            any(grepl(x, pattern="mf[[timetrans$var[i]]]", fixed=TRUE)))) + 1,
          print=FALSE,
          tracer = quote({
            mf <- na.action(mf)
            omit <- attr(mf, "na.action")
            if (!is.null(omit)) {
              Y  <- Y[-omit,]
              strats <- strats[-omit]
            }
          })
  ))
  on.exit({
    suppressMessages(try(untrace(coxph), silent = TRUE))
  })
  
  # Setup call and fit model
  newfrml <- formula(paste(c(newfrml, paste(newtrmstrings, collapse="+"))))
  environment(newfrml) <- newfrmlenv
  pcoxdata <- list2df(as.list(newfrmlenv))
  datameans <- sapply(as.list(newfrmlenv), mean)
  newcall <- expand.call(pcox, call)
  newcall$fitter <- type <- newcall$bs.int <- newcall$bs.yindex <-
    newcall$fitter <- newcall$method <- newcall$eps <- newcall$knots <- NULL
  newcall$formula <- newfrml
  newcall$x <- x
  tt.funcs <- t.funcs[t.types=="tt" & !is.na(t.types)]
  if (length(tt.funcs)) newcall$tt <- quote(tt.funcs)
  newcall$na.action <- na.omit_pcox
  newcall$data <- quote(pcoxdata)
  newcall <- modify_call(newcall, dots)
  newcall$fitter <- newcall$tensortype <- NULL
  newcall[[1]] <- fitter
  res <- eval(newcall)
  
  # Map paramters
  first.para <- 1
  for (i in 1:length(smooth)) {
    if (!is.null(smooth[[i]])) {
      for (j in 1:length(smooth[[i]])) {
        nc <- ncol(smooth[[i]][[j]]$X)
        last.para <- first.para + nc - 1
        names(res$coefficients)[first.para:last.para] <- 
          paste(smooth[[i]][[j]]$label, 1:nc, sep=".")
        names(smooth[[i]][j]) <- smooth[[i]][[j]]$label
        smooth[[i]][[j]]$first.para <- first.para
        smooth[[i]][[j]]$last.para  <- last.para
        first.para <- last.para + 1
      }
      paramap[[i]] <- do.call("c", lapply(smooth[[i]], function(sm)
        sm$first.para:sm$last.para))
    } else {
      paramap[[i]] <- first.para
      first.para <- first.para + 1
    }
  }
  
  # Add smooth objects to the environment of the t.funcs
  sapply(1:length(smooth), function(i) {
    if (!is.null(smooth[[i]])) {
      environment(t.funcs[[i]])$smooth <- smooth[[i]]
    }
    invisible(NULL)
  })
  
  trmmap <- newtrmstrings
  labelmap <- sapply(smooth, function(x) x[[1]]$label)
  names(labelmap) = names(trmmap) <- names(terms)
  # We really don't need trmmap and labelmap? At least not like this.....
  
  # Create smoothmap
  cnt <- 1
  sm.length <- sapply(smooth, length)
  sm.cumsum <- cumsum(sm.length)
  smoothmap <- lapply(1:length(smooth), function(i) {
    if (sm.length[i]>0)
      (sm.cumsum[i]-sm.length[i]+1):(sm.cumsum[i])
  })
  smooth <- do.call("c", smooth)
  smoothdata <- smoothdata[!sapply(smoothdata, is.null)]
  
  # Create varlst
  varlst <- do.call("c", lapply(varmap, function(x) {
    if (is.list(x)) names(x)
    else x
  }))
  
  # Create termtype (type of special term)
  termtype <- rep("par",length(terms))
  for (i in 1:length(specials)) termtype[specials[[i]]-1] <- names(specials)[i]
  
  ret <- list(call=call, formula = formula, method = method,
              responsename = responsename, surv = surv,
              termtype=termtype, termmap = trmmap, labelmap = labelmap,
              varmap = varmap, smoothmap = smoothmap, paramap = paramap,
              t.funcs = t.funcs, t.types = t.types, smooth = smooth,
              smoothdata = smoothdata,
              datameans = datameans, terms=tf)
  res$pcox <- ret
  class(res) <- c("pcox", class(res))
  res
}
