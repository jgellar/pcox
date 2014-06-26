
control.aic <- function (parms, iter, old, n, df, loglik) {
  if (iter == 0) {
    if (is.null(parms$init)) 
      theta <- 0.005
    else theta <- parms$init[1]
    return(list(theta = theta, done = FALSE))
  }
  if (length(parms$type)) 
    type <- parms$type
  else type <- "aic"
  if (n < df + 2) 
    dfc <- (df - n) + (df + 1) * df/2 - 1
  else dfc <- -1 + (df + 1)/(1 - ((df + 2)/n))
  if (iter == 1) {
    history <- c(theta = old$theta, loglik = loglik, df = df, 
                 aic = loglik - df, aicc = loglik - dfc, epic = loglik-2*df)
    if (length(parms$init) < 2) 
      theta <- 1
    else theta <- parms$init[2]
    temp <- list(theta = theta, done = FALSE, history = history)
    return(temp)
  }
  history <- rbind(old$history, c(old$theta, loglik, df, loglik - 
                                    df, loglik - dfc, loglik - 2*df))
  if (is.null(parms$trace)) 
    trace <- FALSE
  else trace <- parms$trace
  if (iter == 2) {
    theta <- mean(history[, 1])
    return(list(theta = theta, done = FALSE, history = history, tst=4))
  }
  
  if (type=="caic") {
    aic <- history[,5]
  } else if (type=="epic") {
    aic <- history[,6]
  } else {
    aic <- history[,4]
  }
  # if (correct) 
  # aic <- history[, 5]
  # else aic <- history[, 4]
  done <- (abs(1 - aic[iter]/aic[iter - 1]) < parms$eps)
  x <- history[, 1]
  if (x[iter] == max(aic) && x[iter] == max(x)) 
    newtheta <- 2 * max(x)
  else newtheta <- survival:::frailty.brent(x, aic, lower = parms$lower, 
                                            upper = parms$upper)
  if (length(parms$trace) && parms$trace) {
    print(history)
    cat("    new theta=", format(newtheta), "\n\n")
  }
  list(theta = newtheta, done = done, history = history, tst=4)
}

