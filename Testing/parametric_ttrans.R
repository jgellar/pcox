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

object <- s(s, t, bs="pb", xt=list(bs="tp", tf=t_transforms))
sm1 <- smooth.construct.pb.smooth.spec(object, data = data, knots=NULL)
sm2 <- smoothCon(object, data=data, knots=NULL, absorb.cons=TRUE)[[1]]

sm3 <- smoothCon(s(s, y, t, bs="pb", xt=list(bs="tp", tf=t_transforms)),
                 data=data, knots=NULL, absorb.cons=TRUE)[[1]]


sm4 <- smoothCon(s(s,t), data=data, knots=NULL, absorb.cons=TRUE)[[1]]

pm1 <- PredictMat(sm1, data=data)
pm2 <- PredictMat(sm2, data=data)
pm3 <- PredictMat(sm3, data=data)



m <- gam(y ~  s(s, t, bs="pb", xt=list(bs="tp", tf=t_transforms)),
  data=data)
m2 <- gam(y ~  s(s, t, bs="pb", xt=t_transforms, k=15),
         data=data)
data$fit <- fitted(m)
data$fit2 <- fitted(m2)
ggplot2::qplot(x=s, y=t, fill = f, data=data, geom="tile")
ggplot2::qplot(x=s, y=t, fill = fit2, data=data, geom="tile")




# Now try w functional variables



