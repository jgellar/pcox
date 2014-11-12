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
#' @import mgcv refund coxme survival
#' @export
#' @author Jonathan Gellar <jgellar1@@jhu.edu>
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
pcox2 <- function(formula, data, method=c("aic","caic","epic","reml", "ml", "fixed", "df"),
                 eps=1e-6, knots=NULL, ...) {
  # Preliminaries...
  call <- match.call()
  method <- match.arg(method)
  dots <- list(...)
  fitter <- ifelse(method %in% c("reml","ml"), "coxme", "coxph")
  if (length(dots)) {
    validDots <- ifelse(fitter=="coxme", names(formals(coxme)),
                        names(formals(coxph)))
    notUsed <- names(dots)[!(names(dots) %in% validDots)]
    if (length(notUsed)) 
      warning("Arguments <", paste(notUsed, collapse = ", "), 
              "> supplied but not used.")
  }  
  
  # Organize terms
  tf <- terms.formula(formula, specials = c("tv", "s", "te", "t2", 
                                            "bf", "hf", "cf",
                                            "strata", "cluster", "frailty"))
  trmstrings <- attr(tf, "term.labels")
  terms <- sapply(trmstrings, function(trm) as.call(parse(text = trm))[[1]], 
                  simplify = FALSE)
  frmlenv <- environment(formula)
  specials <- attr(tf, "specials")
  where.s  <- specials$s - 1
  where.bf <- specials$bf - 1
  where.hf <- specials$hf - 1
  where.cf <- specials$cf - 1
  
  where.sm <- c(where.s, where.bf, where.hf, where.cf)
  where.par <- if (length(trmstrings)) {
    which(!(1:length(trmstrings) %in% where.sm))
  } else {
    numeric(0)
  }
  
  # Organize response and environment
  responsename <- attr(tf, "variables")[2][[1]]
  newfrml <- paste(safeDeparse(responsename), "~", sep = "")
  newfrmlenv <- new.env()
  evalenv <- if ("data" %in% names(call)) 
    eval(call$data)
  else NULL
  
  nobs <- nrow(eval(responsename, envir = evalenv, enclos = frmlenv))
  fitter <- as.symbol(fitter)
  # Need to be modified for survival object?
  if (is.call(responsename)) {
    # responsename is a call to Surv, assign its arguments to newfrmlenv
    sapply(as.list(responsename[-1]), function(x) {
      assign(x = deparse(x),
             value = eval(x, envir=evalenv, enclos=frmlenv),
             envir = newfrmlenv)
      invisible(NULL)
    })
  } else {
    # assign responsename directly to newfrmlenv
    assign(x = deparse(responsename),
           value = eval(responsename, envir = evalenv, enclos = frmlenv),
           envir = newfrmlenv)
  }
  newtrmstrings <- attr(tf, "term.labels")
  tt_funcs <- vector(mode = "list", length = length(trmstrings))
  smooth = specs <- vector("list", length=length(newtrmstrings))
  
  
  #################
  # Process Terms #
  #################
  
  # Smooth terms
  if (length(where.sm)) {
    for (i in where.sm) {
    #sapply(where.sm, function(i) {
      # Evaluate term to get data and tt function
      trm <- eval(terms[[i]], envir=evalenv, enclos=frmlenv)
      
      # Extract tt function, assign env and index, and save
      #tt.i <- attr(trm, "tt")
      tt.i <- trm$tt
      environment(tt.i)$env <- environment() # Or new environment?
      environment(tt.i)$index <- i
      environment(tt.i)$method <- method
      environment(tt.i)$eps <- eps
      #attr(trm, "tt") <- NULL
      tt_funcs[[i]] <- tt.i
      
      # Assign term to newfrmlenv and update newtrmstrings
      nm <- paste0("smooth",i)
      assign(x=nm, trm$x, envir=newfrmlenv)
      newtrmstrings[i] <- paste0("tt(",nm,")")
      
      #invisible(NULL)
    #})
    }
  }
  
  # Parametric terms
  if (length(where.par)) {
    if ("data" %in% names(call))
      frmlenv <- list2env(eval(call$data), frmlenv)
    lapply(terms[where.par], function(x) {
      nms <- if (!is.null(names(x))) {
        all.vars(x[names(x) == ""])
      }
      else all.vars(x)
      sapply(nms, function(nm) {
        stopifnot(length(get(nm, envir = frmlenv)) == nobs)
        assign(x = nm, value = get(nm, envir = frmlenv), envir = newfrmlenv)
        invisible(NULL)
      })
      invisible(NULL)
    })
    tt_funcs[[where.par]] <- function(x, ...) x
  }  
  

  # Setup call and fit model
  newfrml <- formula(paste(c(newfrml, paste(newtrmstrings, collapse="+"))))
  environment(newfrml) <- newfrmlenv
  pcoxdata <- list2df(as.list(newfrmlenv))
  datameans <- sapply(as.list(newfrmlenv), mean)
  newcall <- expand.call(pcox, call)
  newcall$fitter <- type <- newcall$bs.int <- newcall$bs.yindex <- newcall$fitter <-
    newcall$method <- newcall$eps <- newcall$knots <- NULL
  newcall$formula <- newfrml
  newcall$tt <- quote(tt_funcs)
  newcall$data <- quote(pcoxdata)
  newcall$fitter <- newcall$tensortype <- NULL
  newcall[[1]] <- fitter
  res <- eval(newcall)
  
  # Clean up
  if (length(where.sm)) {
    first.para <- 1
    for (i in 1:length(smooth)) {
      if (!is.null(smooth[[i]])) {
        environment(tt_funcs[[i]])$sm <- smooth[[i]]
        nc <- ncol(smooth[[i]]$X)
        last.para <- first.para+nc-1
        names(res$coef)[first.para:last.para] <- paste(smooth[[i]]$label, 1:nc, sep=".")
        names(smooth)[i] <- smooth[[i]]$label
        smooth[[i]]$first.para <- first.para
        smooth[[i]]$last.para  <- last.para
        first.para <- last.para+1
      } else {
        first.para <- first.para + 1
      }
    }
  }
  
  trmmap <- newtrmstrings
  names(trmmap) <- names(terms)
  labelmap <- as.list(trmmap)
  lbls <- sapply(smooth, function(x) x$label)
  
  smooth[sapply(smooth,is.null)] <- NULL
  specs[sapply(specs,is.null)] <- NULL
  
  termtype <- rep("par",length(terms))
  for (i in 1:length(specials)) termtype[specials[[i]]-1] <- names(specials)[i]
  
  ret <- list(formula = formula, method = method, responsename = responsename, nobs = nobs,
              termtype=termtype, termmap = trmmap, labelmap = labelmap,
              where = list(where.s=where.s, where.bf=where.bf, where.hf=where.hf,
                           where.cf=where.cf, where.par = where.par),
              datameans = datameans, specs=specs, smooth=smooth, terms=tf)
  res$pcox <- ret
  class(res) <- c("pcox", class(res))
  res
}

