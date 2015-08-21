#' Fit a Penalized Cox Regression Model
#' 
#' This function fits a penalized Cox regression model using penalized splines 
#' to model the smooth coefficients. It is a wrapper for \code{survival::coxph}
# and \code{coxme::coxme}
#' to fit Cox models with smooth covariate effects for
#' (possibly censored) time-to-event data, of the general form \eqn{log[h_i(t)]
#' = log[h_0(t)] + f_1(z_i, t) + f_2(X_i(s), t)} where \eqn{h_i(t)} is the
#' hazard function for subject i, \eqn{h_0(t)} is the baseline hazard function,
#' \eqn{z_i} are scalar covariates, and \eqn{X_i(s)} are functional covariates.
#' Penalized covariates are set up using functions from the \code{mgcv} package
#' to define basis functions and penalization.
#' 
#' @param formula a formula, with special terms indicating penalized covariates;
#'   see details.
#' @param data an optional data frame containing the variables in the model. If
#'   not found in \code{data}, the variables are taken from
#'   \code{environment(formula)}, typcially he environment from which
#'   \code{pcox} is called.
#' @param method method used to optimize the smoothing parameter(s). \code{"aic"}
#'   is the Akaike Informaion Criterion, \code{"caic"} is the corrected AIC
#'   of Hurvich, et al (1998), and \code{"epic"} is the EPIC criterion of
#'   Shinohara, et al (2011).
#' @param eps convergence level for the criterion indicated by \code{method}.
#' @param knots optional list containing the user-specified knot values for
#'   basis construction, as in \code{\link[mgcv]{gam}}.
#' @param x should the design matrix be included in the return object?
#'   Defaults to \code{TRUE}.
#' @param ...   additional arguments for \code{\link[survival]{coxph}} or
#'   \code{\link[coxme]{coxme}}.
#'   
#' @details
#'  Penalized covariates are specified as special terms in the model formula.
#'  The following terms are allowed:
#'  \enumerate{
#'  \item{\code{\link{p}()}: Specifies a generic penalized term. This is primarily
#'    used for smooth functions of scalar covariates (e.g., \eqn{f(x_i)}) or
#'    time-varying effects of scalar covariates (\eqn{\beta(t) x_i})}. However,
#'    all penalized terms in \code{pcox} formula may be specified using \code{p}
#'    with the appropriate \code{limits} argument. In fact, the other special
#'    terms below are simply wrappers for \code{p}.
#'  \item{\code{\link{cf}()}: Wrapper for \code{p} to specify a concurrent or
#'    lagged effect of a time-varying covariate (concurrent functional
#'    predictor), e.g. \eqn{\beta x_i(t)} or \eqn{\beta x_i(t-\delta_0)}.}
#'  \item{\code{\link{bf}()}: Wrapper for \code{p} to specify a baseline
#'    functional predictor, e.g. \eqn{\int X_i(s)\beta(s) ds}.}
#'  \item{\code{\link{hf}()}: Wrapper for \code{p} to specify a historical
#'    functional predictor, e.g. \eqn{\int_0^t X_i(s)\beta(s,t) ds}.}
#'  }
#'  See the help files for each of the above terms for details regarding usage.
#'  In addition to these special terms, \code{\link{frailty}},
#'  \code{\link{strata}}, \code{\link{cluster}}, and \code{\link{tt}} terms are
#'  allowed, as in a \code{coxph} formula.
#'  
#' @importFrom mgcv gam gam.fit s te t2
#' @importFrom survival coxph Surv
#' @importFrom pryr modify_call
#' @export
#' @author Jonathan Gellar <JGellar@@mathematica-mpr.com> and Fabian Scheipl
#' 
#' @return a fitted \code{pcox} object. This is a \code{coxph} object with
#' additional information in the \code{pcox} entry.
#'
#' @references
#' Gellar, Jonathan E., Elizabeth Colantuoni, Dale M. Needham, and Ciprian M.
#' Crainiceanu (2015). Cox Regression Models with Functional Covariates for Survival
#' Data. Statistical Modeling, 15 (3): 256-278.
#' 
#' Hurvich, C. M., Simonoff, J. S., and Tsai, C. L. (1998). Smoothing Parameter
#' Selection in Nonparametric Regression Using an Improved Akaike Information
#' Criterion. Journal of the Royal Statistical Society. Series B
#' (Methodological): 60(2), 271-293.
#' 
#' Shinohara, R., Crainiceanu, C., Caffo, B., and Reich, D. (2011). Longitudinal
#' Analysis of Spatiotemporal Processes : A Case Study of Dynamic Contrast-Enhanced
#' Magnetic Resonance Imaging in Multiple Sclerosis. Retrieved from
#' http://biostats.bepress.com/jhubiostat/paper231/.
#' 
#' @seealso \code{mgcv}'s \code{\link[mgcv]{smooth.terms}} for details of 
#'   \code{mgcv} syntax and available spline bases and penalties; the related 
#'   \code{\link[refund]{pfr}} and \code{\link[refund]{pffr}} from 
#'   \code{refund}.
#'   
#' @examples
#' # Generate some data
#' N <- 500
#' J <- 200
#' x <- runif(N, 0, 2*pi)
#' male <- rbinom(N, size = 1, prob=.5)
#' 
#' # Simple linear terms (can be done with coxph)
#' eta1 <- matrix(.6*x + .75*male, nrow=N, ncol=J)
#' dat1 <- simTVSurv(eta1, data.frame(x=x, male=male))
#' fit1 <- pcox(Surv(time, event) ~ x + male, data=dat1)
#' fit1b<- survival::coxph(survival::Surv(time,event) ~ x + male, data=dat1)
#' all.equal(fit1$coefficients, fit1b$coefficients)
#' 
#' # Smooth function of scalar x
#' eta2 <- matrix(sin(x) + .75*male, nrow=N, ncol=J)
#' dat2 <- simTVSurv(eta2, data.frame(x=x, male=male))
#' fit2 <- pcox(Surv(time, event) ~ p(x, linear=FALSE) + male, data=dat2)
#' est2 <- coef(fit2)
#' #plot(value ~ x, data=est2, type="l", lwd=1.5, ylim=1.25*range(est2$value))
#' #lines(sin(seq(0,2*pi,by=.1)) ~ seq(0,2*pi,by=.1), col="red")
#' #lines(value + 1.96*se ~ x, data=est2, lty=2)
#' #lines(value - 1.96*se ~ x, data=est2, lty=2)
#' #legend("topright", c("Estimate", "CI", "Truth"),
#' #       lty=c(1,2,1), col=c(1,1,2))
#' 
#' # See help files for p, cf, bf, and hf for more examples

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
  dots$tt <- NULL
  fitter <- ifelse(method %in% c("reml","ml"), "coxme", "coxph")
  
  # Organize terms
  tf <- terms.formula(formula, specials = c("p", "bf", "hf", "cf", "tt"))
  
  #tf <- terms.formula(formula, specials = c("p", "bf", "hf", "cf",
  #                                          "strata", "cluster", "frailty"))
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
  t.funcs = varmap = smooth = smoothdata  <-
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
      nms <- all.vars(term.i[[1]])
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
            tmp <- attr(mf, "na.action")
            attr(mf, "na.action") <- NULL
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
  
  # Finish seting up tt functions
  tt.funcs <- t.funcs[t.types=="tt" & !is.na(t.types)]
  if (length(tt.funcs)) {
    where.tt <- specials$tt - 1
    if (length(where.tt)) {
      # Insert old tt's in correct place
      old.tt <- list(...)$tt
      if (is.function(old.tt))
        old.tt <- list(old.tt)
      tt.locs <- which(t.types=="tt")
      # tt.locs should always be the same length as tt.funcs
      for (i in 1:length(where.tt)) {
        locs <- which(tt.locs < where.tt[i])
        tt.funcs <- append(tt.funcs, old.tt[[i]], 
                           ifelse(length(locs), max(locs), 0))
        t.types[where.tt] <- "tt"
        tt.locs <- which(t.types=="tt") # Adds new tt.location
        #newtt.locs[newtt.locs>where.tt[i]] <-
        #  newtt.locs[newtt.locs>where.tt[i]] + 1
      }
    }
  }
  
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
  if (length(tt.funcs)) newcall$tt <- quote(tt.funcs)
  newcall$na.action <- na.omit_pcox
  newcall$data <- quote(pcoxdata)
  newcall <- modify_call(newcall, dots)
  newcall$fitter <- newcall$tensortype <- NULL
  newcall[[1]] <- fitter
  res <- eval(newcall)
  
  # Map smooths
  for (i in 1:length(smooth)) {
    if (!is.null(smooth[[i]])) {
      idxs <- res$assign2[[i]]
      start <- 1
      for (j in 1:length(smooth[[i]])) {
        idxs.j <- idxs[start:(start+ncol(smooth[[i]][[j]]$X)-1)]
        names(res$coefficients)[idxs.j] <- 
          paste(smooth[[i]][[j]]$label, 1:length(idxs.j), sep=".")
        smooth[[i]][[j]]$first.para <- min(idxs.j)
        smooth[[i]][[j]]$last.para  <- max(idxs.j)
        start <- start + length(idxs.j)
      }
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
              varmap = varmap, smoothmap = smoothmap,
              t.funcs = t.funcs, t.types = t.types, smooth = smooth,
              smoothdata = smoothdata,
              datameans = datameans, terms=tf)
  res$pcox <- ret
  class(res) <- c("pcox", class(res))
  res
}
