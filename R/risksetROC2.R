#' risksetROC 2nd version
#' 
#' @keywords internal
#' @importFrom risksetROC CoxWeights
#' @importFrom graphics abline
risksetROC2 <- function (Stime, entry = NULL, status, marker, predict.time, 
          span = NULL, order = 1, window = "asymmetric", 
          prop = 0.5, plot = TRUE, type = "l", xlab = "FP", ylab = "TP", 
          ...) 
{
  p = 1
  eta = marker
  if (length(entry) == 0) {
    entry = rep(0, NROW(Stime))
  }
  time = entry
  time2 = Stime
  new.eta <- eta
  
  out = CoxWeights(marker = new.eta, entry = time, status = status, 
                   predict.time = predict.time, Stime = time2)
  if (plot == TRUE) {
    plot(out$FP, out$TP, type = type, xlab = xlab, ylab = ylab, 
         ...)
    abline(c(0, 0), c(1, 1))
  }
  return(out)
}
