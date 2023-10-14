
## probably need to remove sigquant and sigdf at some point.

RMST_BART <- function(
    U, delta, x.train, Gweights, x.test=matrix(0.0,0,0),
    tau=NULL, sgrid=NULL, sigma.mu=NULL,
    transformation="identity",
    sparse=FALSE, theta=0, omega=1,
    a=0.5, b=1, augment=FALSE, rho=NULL,
    xinfo=matrix(0.0,0,0), usequants=TRUE,
    cont=FALSE, rm.const=TRUE,
    k=2.0, power=2.0, base=.95, sigmaf=NA,
    ntree=200L, numcut=100L,
    ndpost=1000L, nskip=100L, keepevery=1L,
    nkeeptrain=ndpost, nkeeptest=ndpost,
    nkeeptestmean=ndpost, nkeeptreedraws=ndpost,
    printevery=100L) {
  ## xmat is the matrix that only contains the rows of
  ## X where delta=1
  if(ncol(x.train) > 1) {
    xmat <- x.train[delta==1,]
  } else{
    xmat <- matrix(x.train[delta==1,], nrow=sum(delta==1), ncol=1)
    colnames(xmat) <- colnames(X)
  }

  kk <- k
  temp = bartModelMatrix(xmat, numcut, usequants=usequants,
                         cont=cont, xinfo=xinfo, rm.const=rm.const)
  x.train = t(temp$X)
  numcut = temp$numcut
  xinfo = temp$xinfo
  if(length(x.test)>0) {
    x.test = bartModelMatrix(x.test)
    x.test = t(x.test[ , temp$rm.const])
  }
  rm.const <- temp$rm.const
  grp <- temp$grp
  rm(temp)

  n <- ncol(x.train)
  p <- nrow(x.train)
  np <- ncol(x.test)
  if(length(rho)==0) rho=p
  if(length(rm.const)==0) rm.const <- 1:p
  if(length(grp)==0) grp <- 1:p

  Gweights <- c(t(Gweights))
  if(length(Gweights) != n*(ndpost*keepevery + nskip + 1)) {
      stop("Gweights does not have the correct dimensions")
  }
  ## As a default, set tau to maximum of observed survival times
  if(is.null(tau)) {
    ## what to do if all delta==0?
    tau <- max(U[delta==1])
  }
  U_tau <- pmin(U[delta==1], tau)

  if(is.null(sgrid)) {
    sgrid <- seq(0, tau, length.out=100)
  }

  ## Get KM estimate of censoring distribution and KM inverse censoring weights
  KM_cens <- survfit(Surv(U, 1 - delta) ~ 1)
  GKMfn <- stepfun(c(0, KM_cens$time), c(1, KM_cens$surv, min(KM_cens$surv)))
  GKM_weights <- GKMfn(U_tau)
  if(sum(GKM_weights == 0) > 0) {
     GKM_weights[GKM_weights==0] <- min(KM_cens$surv[KM_cens$surv > 0])
  }

  ## Define centered version of U_tau
  ## and set sigma.mu if its values was not input to the function
  if(transformation=="identity") {
    ## compute muhatb
    muhatb <- mean(U_tau/GKM_weights)
    Y_tau <- U_tau - muhatb
    if(is.null(sigma.mu)) {
      Ymin <- min(U_tau)
      sigma.mu <- (tau - muhatb - Ymin)/(2*kk*sqrt(ntree))
    }
  } else if(transformation=="log") {
    ## compute muhatb
    muhatb <- mean(log(U_tau)/GKM_weights)
    Y_tau <- log(U_tau) - muhatb
    if(is.null(sigma.mu)) {
      Ymin <- min(U_tau)
      sigma.mu <- (log(tau) - muhatb - log(Ymin))/(2*kk*sqrt(ntree))
    }
  }


  ptm <- proc.time()
  #call
  ## is a vector of weights here.
  nu <- 3
  lambda <- 1.0
  sigest <- 1.0

  res = .Call("cwbart",
              n,  #number of observations in training data
              p,  #dimension of x
              np, #number of observations in test data
              x.train,   #pxn training data x
              Y_tau,
              x.test,   #p*np test data x
              ntree,
              numcut,
              ndpost*keepevery, #this is nd
              nskip, # this is burn
              power,
              base,
              sigma.mu, # should change variable name in C++ file.
              nu,
              lambda,
              sigest,
              Gweights,
              sparse,
              theta,
              omega,
              grp,
              a,
              b,
              rho,
              augment,
              nkeeptrain,
              nkeeptest,
              nkeeptestmean,
              nkeeptreedraws,
              printevery,
              xinfo
  )

  res$mu = muhatb

  if(transformation=="identity") {
     res$yhat.train <- pmax(pmin(res$yhat.train + muhatb, tau), 0.0)
     res$yhat.test <- pmax(pmin(res$yhat.test + muhatb, tau), 0.0)
     res$yhat.train.mean <- colMeans(res$yhat.train)
     res$yhat.test.mean <- colMeans(res$yhat.test)
  } else if(transformation=="log") {
     res$yhat.train <- pmin(res$yhat.train + muhatb, log(tau))
     res$yhat.test <- pmin(res$yhat.test + muhatb, log(tau))
     res$yhat.train.mean <- colMeans(res$yhat.train)
     res$yhat.test.mean <- colMeans(res$yhat.test)
  }
  if(nkeeptreedraws>0)
    names(res$treedraws$cutpoints) = dimnames(x.train)[[1]]
  dimnames(res$varcount)[[2]] = as.list(dimnames(x.train)[[1]])
  dimnames(res$varprob)[[2]] = as.list(dimnames(x.train)[[1]])
  res$varcount.mean <- apply(res$varcount, 2, mean)
  res$varprob.mean <- apply(res$varprob, 2, mean)
  res$rm.const <- rm.const
  attr(res, 'class') <- 'bart_rmst'
  return(res)
}
