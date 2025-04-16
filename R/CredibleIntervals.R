CredibleIntervals <- function(obj, level=0.95) {
  if(level >= 1.0 | level <= 0.0) {
       stop("level must be between 0 and 1")
  }
  alpha_value <- (1 - level)/2
  CredInterval <- matrix(NA, nrow=ncol(obj$yhat.train), ncol=2)
  CredInterval[,1] <- apply(obj$yhat.train, 2, function(x) quantile(x, probs=alpha_value))
  CredInterval[,2] <- apply(obj$yhat.train, 2, function(x) quantile(x, probs=1 - alpha_value))
  return(CredInterval)
}
