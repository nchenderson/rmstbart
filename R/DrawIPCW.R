LambdaCDF <- function(t, nevent, nrisk, c0, delta_alpha) {

  tmp <- (((-1)^(0:nevent))*choose(nevent, 0:nevent))/((nrisk - nevent + c0 + 0:nevent)^(c0*delta_alpha))
  ww <- tmp/sum(tmp)
  ff <- pgamma(t, shape=c0*delta_alpha, rate=nrisk - nevent + c0 + 0:nevent)
  FF <- sum(ww*pgamma(t, shape=c0*delta_alpha, rate=nrisk - nevent + c0 + 0:nevent))
  return(FF)
}

InvLambdaCDF <- function(u, nevent, nrisk, c0, delta_alpha) {
  ff <- function(tt, u, nevent, nrisk, c0, delta_alpha) {
    return( LambdaCDF(tt, nevent, nrisk, c0, delta_alpha) - u)
  }
  ## need a better way to figure out ulimit
  ulimit <- 100
  print(c(nrisk, nevent))
  print(LambdaCDF(0, nevent, nrisk, c0, delta_alpha))
  print(LambdaCDF(100, nevent, nrisk, c0, delta_alpha))

  ans <- uniroot(ff, lower=0, upper=ulimit,
                 u = u, nevent=nevent, nrisk=nrisk, c0=c0, delta_alpha=delta_alpha)$root
  return(ans)
}

DrawLambdas <- function(nevents, nrisks, kappa0, delta_alpha) {
  nbins <- length(nevents)
  lambda_draw <- rep(NA, nbins)
  eta <- kappa0*delta_alpha
  bbeta <- nrisks - nevents + kappa0
  for(k in 1:nbins) {

    if(eta==1) {
      rax <- rbeta(1, shape1=bbeta[k], shape2=nevents[k] + 1)
      lambda_draw[k] <- -log(rax)
    } else {
      done <- FALSE
      while(!done) {
        theta <- rgamma(1, shape=eta + nevents[k], rate=bbeta[k])
        u <- runif(1)
        thresh <- log1p(-exp(-theta)) - log(theta)
        #print(c(thresh, nevents[k]*thresh, log(u)))
        # This approach doesn't work too well if
        #  nevents[k] is large
        if(log(u) < nevents[k]*thresh) {
          lambda_draw[k] <- theta
          done <- TRUE
        }
      }
    }
  }
  return(lambda_draw)
}

DrawIPCW <- function(U, delta, Utau, sgrid, kappa0, delta_alpha) {
  J <- length(sgrid)
  E <- R <- rep(NA, J-1)
  for(j in 2:J) {
    E[j-1] <- sum((1 - delta)*(U > sgrid[j-1])*(U <= sgrid[j]))
    R[j-1] <- sum(U > sgrid[j-1])
  }
  lambda.draw <- DrawLambdas(nevents=E, nrisks=R, kappa0=kappa0, delta_alpha=delta_alpha)
  CumHazFn <- approxfun(sgrid, cumsum(c(0, lambda.draw)))
  weight_ans <- exp(CumHazFn(Utau))
  return(weight_ans)
}


## Testing it out
#n <- 500
#TT <- rexp(n)
#CC <- rexp(n, rate=2.5)
#sgrid <- seq(0, 10, length.out=500)

#U <- pmin(TT, CC)
#delta <- ifelse(TT <= CC, 1, 0)

#tt <- DrawIPCW(U=U, delta=delta, Utau=U, sgrid=sgrid,
#               kappa0=1, delta_alpha=1)

#xx <- seq(0, 3, length.out=1000)
#plot(U[order(tt)], 1/tt[order(tt)])
#lines(xx, exp(-2.5*xx), lwd=3, col="red")

#kappa0 <- 1
#delta_alpha <- 1/kappa0
#Gmat <- matrix(1, nrow=ndraws + burnIn + 1, ncol=length(U_tau))
#for(k in 1:(ndraws + burnIn + 1)) {
#  Gmat[k,] <- DrawIPCW(U=Y.train, delta=delta.train, Utau=U_tau, sgrid=sgrid,
#                       kappa0=kappa0, delta_alpha=delta_alpha)
#}
#Gmat_orig <- 1/sqrt(Gmat)
#Gmeans <- colMeans(1/Gmat)
#Gmat <- sqrt(2*eta_hat_vals[u])*Gmat_orig


