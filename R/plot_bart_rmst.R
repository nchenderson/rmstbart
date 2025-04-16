plot.bart_rmst <- function(x, xlab=NULL, main=NULL, level=0.95) {

   if(level >= 1.0 | level <= 0.0) {
      stop("level must be between 0 and 1")
   }
   alpha_value <- (1 - level)/2
   if(is.null(xlab)) {
       xlab <- "RMST"
   }
   SummaryStat <- matrix(NA, nrow=length(x$yhat.train.mean), ncol=3)
   SummaryStat[,1] <- x$yhat.train.mean
   SummaryStat[,2] <- apply(x$yhat.train, 2, function(x) quantile(x, probs=alpha_value))
   SummaryStat[,3] <- apply(x$yhat.train, 2, function(x) quantile(x, probs=1 - alpha_value))
   SummaryStat <- SummaryStat[order(x$yhat.train.mean),]
   nn <- nrow(SummaryStat)
   xlim_low <- max(min(SummaryStat[,2]), 0)
   xlim_high <- max(SummaryStat[,3])
   plot(0,0, type="n", xlim=c(xlim_low, xlim_high),
        ylim=c(0, nn), xlab=xlab, ylab="Observations sorted by RMST",
        main=main, las=1)
   for(k in 1:nn) {
       points(SummaryStat[k,1],k, pch=16, col="blue")
       lines(c(SummaryStat[k,2], SummaryStat[k,3]), c(k, k), pch=16, col="grey")
   }
}
