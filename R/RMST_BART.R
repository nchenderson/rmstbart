
rmstbart <- function(
    times, delta, x.train, x.test=matrix(0.0,0,0),
    tau=NULL, sgrid=NULL,
    censoring="dependent",
    transformation="identity",
    sparse=FALSE, theta=0, omega=1,
    a=0.5, b=1, augment=FALSE, rho=NULL,
    xinfo=matrix(0.0,0,0), usequants=TRUE, cross.val=TRUE,
    kappa=2.0, power=2.0, base=.95, sigmaf=NA,
    ntree=200L, numcut=100L, ndpost=1000L, nskip=100L,
    keepevery=1L, nkeeptrain=ndpost, nkeeptest=ndpost,
    nkeeptestmean=ndpost, nkeeptreedraws=ndpost,
    printevery=100L) {


    ## This drops constant variables from the analysis
    rm.const <- TRUE
    ## The variable cont states whether or not to assume all
    ## covariates are continuous.
    cont <- FALSE

    ### First check if all observations are censored.
    nd0 <- sum(delta==0)
    if(nd0==nrow(x.train)) {
        stop("All observations are censored")
    }

    ## Save these terms, to use in the censoring model
    x.train.orig <- x.train
    y.train.orig <- times
    delta.orig <- delta
   # cens_dist <- survfit(Surv(y.train.orig, 1-delta.orig) ~ 1)
    #GKM <- stepfun(c(0, cens_dist$time), c(1, cens_dist$surv, min(cens_dist$surv)))


    ## Create matrix xmat. This will be the matrix that only
    ## contains the rows of x.train where events occur (i.e. delta=1)
    if(ncol(x.train) > 1) {
        xmat <- x.train[delta==1,]
        xmat_d0 <- x.train[delta==0,]
    } else{
        xmat <- matrix(x.train[delta==1,], nrow=sum(delta==1), ncol=1)
        colnames(xmat) <- colnames(x.train)

        xmat_d0 <- matrix(x.train[delta==0,], nrow=sum(delta==0), ncol=1)
        colnames(xmat_d0) <- colnames(x.train)
    }

    ## Create model matrix tailored for use by BART:
    temp = bartModelMatrix(xmat, numcut, usequants=usequants,
                           cont=cont, xinfo=xinfo, rm.const=rm.const)

    x.train = t(temp$X)
    numcut = temp$numcut
    xinfo = temp$xinfo
    testset_used <- FALSE
    if(length(x.test)>0 & nd0 > 0) {
        ## If censoring is present, add censored obs to the input
        ## test set
        x.tmp <- rbind(xmat_d0, x.test)
        x.test <- bartModelMatrix(x.tmp)
        x.test <- t(x.test[ , temp$rm.const])
        testset_used <- TRUE
     } else if(length(x.test)==0 & nd0 > 0) {
        ## If censoring is present and no test set is entered,
        ## only treat censoring observations as test set
        x.test <- bartModelMatrix(xmat_d0)
        x.test <- t(x.test[ , temp$rm.const])
     } else if(length(x.test)>0 & nd0 == 0) {
        x.test <- bartModelMatrix(x.test)
        x.test <- t(x.test[ , temp$rm.const])
        testset_used <- TRUE
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
     xmat <- xmat[,rm.const]

     ## As a default, set tau to maximum of observed survival times
     U <- times
     if(is.null(tau)) {
       ## what to do if all delta==0?
       tau <- max(U[delta==1])
     }
     U_tau <- pmin(U[delta==1], tau)

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

        Ymin <- min(U_tau)
        sigma.mu <- (tau - muhatb - Ymin)/(2*kappa*sqrt(ntree))

        ## Fit preliminary AFT model to get candidate tuning parameters
        ## for AFT model.
        AFT_try <- survreg(Surv(exp(y.train.orig), delta.orig) ~ x.train.orig)
        eta_hat <- AFT_try$scale*AFT_try$scale

     } else if(transformation=="log") {
        ## compute muhatb
        muhatb <- mean(log(U_tau)/GKM_weights)
        Y_tau <- log(U_tau) - muhatb

        Ymin <- min(U_tau)
        sigma.mu <- (log(tau) - muhatb - log(Ymin))/(2*kappa*sqrt(ntree))

        ## Fit preliminary AFT model to get candidate tuning parameters
        ## for AFT model.
        AFT_try <- survreg(Surv(y.train.orig, delta.orig) ~ x.train.orig)
        eta_hat <- AFT_try$scale*AFT_try$scale
     }

     cat("Generating censoring weights. \n")
     if(censoring=="dependent") {
          cens_bart <- AFTrees::AFTrees(x.train=x.train.orig, y.train=y.train.orig, status=1-delta.orig,
                                        ndpost=ndpost + nskip, verbose=FALSE)
          Mucens_draws <- cens_bart$m.train
          GmatDep <- matrix(1, nrow=ndpost + nskip + 1, ncol=length(U_tau))
          for(k in 1:(ndpost + nskip)) {
              for(h in 1:length(U_tau)) {
                 log.time.points <- log(U_tau[h])
                 AA <- (log.time.points - cens_bart$locations[k,] - Mucens_draws[k,h])/cens_bart$sigma[k]
                 Cprob <- sum(pnorm(AA, lower.tail=FALSE)*cens_bart$mix.prop[k,])
                 GmatDep[k,h] <- 1/Cprob
              }
           }
           GmatDeporig <- 1/sqrt(GmatDep)
     } else if(censoring=="independent") {
         kappa0 <- 1
         delta_alpha <- 1/kappa0
         Gmat <- matrix(1, nrow=ndpost + nskip + 1, ncol=length(U_tau))
         for(k in 1:(ndraws + burnIn + 1)) {
           Gmat[k,] <- DrawIPCW(U=Y.train, delta=delta.train, Utau=U_tau, sgrid=sgrid,
                                kappa0=kappa0, delta_alpha=delta_alpha)
         }
         GmatDeporig <- 1/sqrt(Gmat)
     }

    # if(length(Gweights) != n*(ndpost*keepevery + nskip + 1)) {
     #  stop("Gweights does not have the correct dimensions")
     #}


     ptm <- proc.time()
     #call
     ## is a vector of weights here.
     nu <- eta_hat
     lambda <- 1.0
     sigest <- 1.0

     nfolds <- 5
     best_eta <- eta_hat
     if(cross.val) {
         eta_hat_vals <- c(0.1*eta_hat, 0.25*eta_hat, 0.5*eta_hat, 0.75*eta_hat, eta_hat)
         ncv <- length(eta_hat_vals)
         cv_scores <- matrix(NA, nrow=ncv, ncol=nfolds)

         if(n < 100) {
            warning("Cross-validation is not recommended if the number of observed events is less than 100.")
         }
     }

     cat("Now, running MCMC for RMST regression function. \n")
     folds <- sample(1:nfolds, size=n, replace=TRUE)
     if(cross.val) {
          for(u in 1:ncv) {
              ## Split into test and train here
              ## use: x.train.tmp, Y_tau.tmp, x.test.tmp
              #GmatDep <- sqrt(2*eta_hat_vals[u])*GmatDeporig
              #Gweights <- c(t(GmatDep))
              for(k in 1:nfolds) {
                  cat("Cross-validation: Fold ", k,"(out of ", nfolds," folds) of setting", u," (out of ", ncv, "settings). \n")
                  Y.train.tmp <- Y_tau[folds!=k]
                  xmat.tmp <- xmat[folds!=k,]
                  xmat.test.tmp <- xmat[folds==k,]
                  Y.test.tmp <- Y_tau[folds==k]

                  temp <- bartModelMatrix(xmat.tmp)

                  x.train.tmp = t(temp)

                  temp <- bartModelMatrix(xmat.test.tmp)
                  x.test.tmp <- t(temp)

                  GmatDep.train <- sqrt(2*eta_hat_vals[u])*GmatDeporig[,folds!=k]
                  Gweights <- c(t(GmatDep.train))
                  ww_dep <- colMeans(1/(GmatDeporig[,folds==k]^2))

                  n.tmp <- ncol(x.train.tmp)
                  np.tmp <- ncol(x.test.tmp)
                  tmp_res <- .Call("cwbart",
                                    n.tmp,  #number of observations in training data
                                    p,  #dimension of x
                                    np.tmp, #number of observations in test data
                                    x.train.tmp,   #pxn training data x
                                    Y.train.tmp,
                                    x.test.tmp,   #p*np test data x
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

                  yhat.tmp <- tmp_res$yhat.test.mean
                  cv_scores[u,k] <- mean(ww_dep*((Y.test.tmp - yhat.tmp)*(Y.test.tmp - yhat.tmp)))
              }
          }
          best_eta_ind <- which.min(rowMeans(cv_scores))
          best_eta <- eta_hat_vals[best_eta_ind]
          ## pick the best eta and run again.
     }
     if(!cross.val) {
         best_eta <- 0.5
         cv_scores <- NULL
     }
     GmatDep <- sqrt(2*best_eta)*GmatDeporig
     Gweights <- c(t(GmatDep))
     cat("Final MCMC run \n")
     tmp_res <- .Call("cwbart",
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

  res <- tmp_res
  res$mu <- muhatb
  if(transformation=="identity") {
     if(testset_used & nd0 > 0) {
         nt <- ncol(res$yhat.test)
         train_draws <- res$yhat.train
         res$yhat.train <- matrix(NA, nrow=nrow(res$yhat.train),
                                  ncol=ncol(res$yhat.train) + nd0)
         res$yhat.train[,delta==1] <- train_draws
         res$yhat.train[,delta==0] <- res$yhat.test[,1:nd0]

         res$yhat.test <- res$yhat.test[,(nd0+1):nt]

         res$yhat.train <- pmax(pmin(res$yhat.train + muhatb, tau), 0.0)
         res$yhat.test <- pmax(pmin(res$yhat.test + muhatb, tau), 0.0)

         res$yhat.train.mean <- colMeans(res$yhat.train)
         res$yhat.test.mean <- colMeans(res$yhat.test)
     } else if(!testset_used & nd0 > 0) {
         train_draws <- res$yhat.train
         res$yhat.train <- matrix(NA, nrow=nrow(res$yhat.train),
                                  ncol=ncol(res$yhat.train) + nd0)
         res$yhat.train[,delta==1] <- train_draws
         res$yhat.train[,delta==0] <- res$yhat.test[,1:nd0]

         res$yhat.train <- pmax(pmin(res$yhat.train + muhatb, tau), 0.0)
         res$yhat.train.mean <- colMeans(res$yhat.train)

         res$yhat.test <- NULL
         res$yhat.test.mean <- NULL
     } else if(testset_used & nd0 == 0) {
         res$yhat.train <- pmax(pmin(res$yhat.train + muhatb, tau), 0.0)
         res$yhat.train.mean <- colMeans(res$yhat.train)

         res$yhat.test <- pmax(pmin(res$yhat.test + muhatb, tau), 0.0)
         res$yhat.test.mean <- colMeans(res$yhat.test)
     } else {
         res$yhat.train <- pmax(pmin(res$yhat.train + muhatb, tau), 0.0)
         res$yhat.train.mean <- colMeans(res$yhat.train)
     }

  } else if(transformation=="log") {
    if(testset_used & nd0 > 0) {
       nt <- ncol(res$yhat.test)
       train_draws <- res$yhat.train
       res$yhat.train <- matrix(NA, nrow=nrow(res$yhat.train),
                               ncol=ncol(res$yhat.train) + nd0)
       res$yhat.train[,delta==1] <- train_draws
       res$yhat.train[,delta==0] <- res$yhat.test[,1:nd0]

       res$yhat.test <- res$yhat.test[,(nd0+1):nt]

       res$yhat.train <- pmin(res$yhat.train + muhatb, log(tau))
       res$yhat.test <- pmin(res$yhat.test + muhatb, log(tau))

       res$yhat.train.mean <- colMeans(res$yhat.train)
       res$yhat.test.mean <- colMeans(res$yhat.test)
    } else if(!testset_used & nd0 > 0) {
       train_draws <- res$yhat.train
       res$yhat.train <- matrix(NA, nrow=nrow(res$yhat.train),
                               ncol=ncol(res$yhat.train) + nd0)
       res$yhat.train[,delta==1] <- train_draws
       res$yhat.train[,delta==0] <- res$yhat.test

       res$yhat.train <- pmin(res$yhat.train + muhatb, log(tau))
       res$yhat.train.mean <- colMeans(res$yhat.train)

       res$yhat.test <- NULL
       res$yhat.test.mean <- NULL
    } else if(testset_used & nd0 == 0) {
       res$yhat.train <- pmin(res$yhat.train + muhatb, log(tau))
       res$yhat.train.mean <- colMeans(res$yhat.train)

       res$yhat.test <- pmin(res$yhat.test + muhatb, tau)
       res$yhat.test.mean <- colMeans(res$yhat.test)
    } else {
       res$yhat.train <- pmin(res$yhat.train + muhatb, log(tau))
       res$yhat.train.mean <- colMeans(res$yhat.train)
    }
  }
  if(nkeeptreedraws>0)
    names(res$treedraws$cutpoints) = dimnames(x.train)[[1]]

  dimnames(res$varcount)[[2]] = as.list(dimnames(x.train)[[1]])
  dimnames(res$varprob)[[2]] = as.list(dimnames(x.train)[[1]])
  res$varcount.mean <- apply(res$varcount, 2, mean)
  res$varprob.mean <- apply(res$varprob, 2, mean)
  res$rm.const <- rm.const
  res$eta <- eta_hat
  res$cv.scores <- cv_scores
  print(dim(GmatDeporig))
  nrG <- nrow(GmatDeporig)
  res$censoring.weights <- GmatDeporig[(nskip + 1):nrG,]
  print(dim(res$censoring.weights))
  attr(res, 'class') <- 'bart_rmst'
  return(res)
}
