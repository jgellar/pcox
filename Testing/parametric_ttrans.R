# generate data y(s,t) = sum_k g_k(t) f_k(s)
library("mgcv")

t_grid <- seq(0, 1, l=50)
s_grid <- seq(0, 1, l=50)
data <- expand.grid(s=s_grid, t=t_grid)

t_transforms <- list(function(t) t^0, function(t) t, function(t) sin(6*pi*t))
s_shapes <- list(function(s) sin(4 * pi * s), function(s) cos(2*pi*s),
  function(s) log1p(s))
data <- within(data, {
  ft <- vapply(t_transforms, function(f) do.call(f, list(t=t)), numeric(length(t)))
  fs <- vapply(s_shapes, function(f) do.call(f, list(s=s)), numeric(length(t)))
  f <- rowSums(ft * fs)
  y <- f + .1 * rnorm(length(f))
})
ggplot2::qplot(x=s, y=t, fill = rowSums(fs), data=data, geom="tile")
ggplot2::qplot(x=s, y=t, fill = rowSums(ft), data=data, geom="tile")
ggplot2::qplot(x=s, y=t, fill = f, data=data, geom="tile")
ggplot2::qplot(x=s, y=t, fill = y, data=data, geom="tile")

#'
#'
object <- s(s, t, bs="parametrictime", xt=list(bs="tp", tt=t_transforms))
smooth.construct.parametrictime.smooth.spec <- function(object,data,knots){
  ## todo: input checks

  #set up design for t
  t_vec <- data[[object$term[2]]]
  n_transform <- length(object$xt$tt)
  n <- length(data[[object$term[2]]])
  tt_name <- paste0(object$term[2], "_tt_", 1:n_transform)
  t_X <-  vapply(object$xt$tt,
    function(f) do.call(f, list(t=data[[object$term[2]]])), numeric(n))

  #set up smooth over s given args in object$xt
  s_args <- object$xt[!grepl("tt", names(object$xt))]
  s_smoothspec <- do.call(mgcv::s, append(as.name(object$term[1]), s_args))
  sm <- smooth.construct(s_smoothspec, data = data, knots = NULL)
  # modify smooth term:

  sm$term <- object$term
  sm$bs.dim <-  sm$bs.dim * n_transform
  sm$null.space.dim <-  sm$null.space.dim * n_transform
  sm$df <- sm$df * n_transform
  sm$dim <- 2
  sm$label <- paste0("f(", object$term[2], ")*", sm$label)
  sm$xt <- object$xt
  sm$X <- mgcv::tensor.prod.model.matrix(list(t_X, sm$X))
  sm$S[[1]] <- diag(n_transform) %x% sm$S[[1]]
  ## ctrl-c-v from smooth.construct.tp.smooth.spec:
  if (sm$drop.null > 0) {
    ind <- 1:(sm$bs.dim - sm$null.space.dim)
    if (FALSE) { ## nat param version
     np <- nat.param(sm$X, sm$S[[1]], rank=sm$bs.dim - sm$null.space.dim, type=0)
     sm$P <- np$P
     sm$S[[1]] <- diag(np$D)
     sm$X <- np$X[,ind]
    } else { ## original param
     sm$S[[1]] <-sm$S[[1]][ind,ind]
     sm$X <-sm$X[, ind]
     sm$cmX <- colMeans(object$X)
     sm$X <- sweep(object$X, 2, object$cmX)
    }
    sm$null.space.dim <- 0
    sm$df <- sm$df - M
    sm$bs.dim <-sm$bs.dim -M
    sm$C <- matrix(0,0,ncol(object$X)) # null constraint matrix
  }
  sm
}

m <- gam(y ~  s(s, t, bs="parametrictime", xt=list(bs="tp", tt=t_transforms)),
  data=data)
data$fit <- fitted(m)
ggplot2::qplot(x=s, y=t, fill = f, data=data, geom="tile")
ggplot2::qplot(x=s, y=t, fill = fit, data=data, geom="tile")
