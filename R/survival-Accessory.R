# These functions are taken directly from the survival package, with
# permission from Terry Therneau. Copied @ 2015-07-06 from survival 2.38-2.

frailty.brent <- function (x, y, lower, upper) {
  n <- length(x)
  if (length(y) != n) 
    stop("Length mismatch for x and y")
  if (n < 3) 
    return(mean(x))
  ord <- order(x)
  xx <- x[ord]
  yy <- y[ord]
  best <- (1:n)[yy == max(y)]
  if (length(best) > 1) 
    stop("Ties for max(y), I surrender")
  if (best == 1) {
    new <- xx[1] - 3 * (xx[2] - xx[1])
    if (!missing(lower) && !is.null(lower) && new < lower) 
      new <- lower + (min(xx[xx > lower]) - lower)/10
    return(new)
  }
  if (best == n) {
    new <- xx[n] + 3 * (xx[n] - xx[n - 1])
    if (!missing(upper) && !is.null(upper) && new > upper) 
      new <- upper + (max(xx[xx < upper]) - upper)/10
    return(new)
  }
  xx <- xx[(best - 1):(best + 1)]
  yy <- yy[(best - 1):(best + 1)]
  temp1 <- (xx[2] - xx[1])^2 * (yy[2] - yy[3]) - (xx[2] - xx[3])^2 * 
    (yy[2] - yy[1])
  temp2 <- (xx[2] - xx[1]) * (yy[2] - yy[3]) - (xx[2] - xx[3]) * 
    (yy[2] - yy[1])
  new <- xx[2] - 0.5 * temp1/temp2
  if (new < xx[1] || new > xx[3] || ((n > 4) && (new - x[n]) > 
                                       0.5 * abs(x[n - 1] - x[n - 2]))) {
    if ((xx[2] - xx[1]) > (xx[3] - xx[2])) 
      return(xx[2] - 0.38 * (xx[2] - xx[1]))
    else return(xx[2] + 0.32 * (xx[3] - xx[2]))
  }
  else return(new)
}

frailty.controldf <- function (parms, iter, old, df) {
  if (iter == 0) {
    theta <- parms$guess
    return(list(theta = theta, done = FALSE, history = cbind(thetas = parms$thetas, 
                                                             dfs = parms$dfs)))
  }
  eps <- parms$eps
  if (length(eps) == 0) 
    eps <- 0.1
  thetas <- c(old$history[, 1], old$theta)
  dfs <- c(old$history[, 2], df)
  nx <- length(thetas)
  if (nx == 2) {
    theta <- thetas[1] + (thetas[2] - thetas[1]) * (parms$df - 
                                                      dfs[1])/(dfs[2] - dfs[1])
    if (parms$df > df) 
      theta <- theta * 1.5
    return(list(theta = theta, done = FALSE, history = cbind(thetas = thetas, 
                                                             dfs = dfs), half = 0))
  }
  else {
    done <- (iter > 1 && (abs(dfs[nx] - parms$df) < eps))
    x <- thetas
    y <- dfs
    target <- parms$df
    if (abs((y[nx] - target)/(y[nx - 1] - target)) > 0.6) 
      doing.well <- FALSE
    else doing.well <- TRUE
    ord <- order(x)
    if ((x[1] - x[2]) * (y[1] - y[2]) > 0) 
      y <- y[ord]
    else {
      y <- -1 * y[ord]
      target <- -target
    }
    x <- x[ord]
    if (all(y > target)) 
      b1 <- 1
    else if (all(y < target)) 
      b1 <- nx - 2
    else {
      b1 <- max((1:nx)[y <= target])
      if (!doing.well && (is.null(old$half) || old$half < 
                            2)) {
        if (length(parms$trace) && parms$trace) {
          print(cbind(thetas = thetas, dfs = dfs))
          cat("  bisect:new theta=", format(mean(x[b1 + 
                                                     0:1])), "\n\n")
        }
        return(list(theta = mean(x[b1 + 0:1]), done = done, 
                    history = cbind(thetas = thetas, dfs = dfs), 
                    half = max(old$half, 0) + 1))
      }
      if ((b1 + 1) == nx || (b1 > 1 && ((target - y[b1]) < 
                                          (y[b1 + 1] - target)))) 
        b1 <- b1 - 1
    }
    b2 <- b1 + 1:2
    xx <- log(x[b2] - x[b1])
    yy <- log(y[b2] - y[b1])
    power <- diff(yy)/diff(xx)
    a <- yy[1] - power * xx[1]
    newx <- (log(target - y[b1]) - a)/power
    if (length(parms$trace) && parms$trace) {
      print(cbind(thetas = thetas, dfs = dfs))
      cat("  new theta=", format(x[b1] + exp(newx)), "\n\n")
    }
    list(theta = x[b1] + exp(newx), done = done, history = cbind(thetas = thetas, 
                                                                 dfs = dfs), half = 0)
  }
}