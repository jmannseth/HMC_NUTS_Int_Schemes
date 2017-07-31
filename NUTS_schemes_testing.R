## Copyright (C) 2015 Cole Monnahan
## License: GPL-2

#' [BETA VERSION] Draw samples from the posterior of a TMB model using a
#' specified MCMC algorithm.
mcmc2 <- function(obj, nsim, algorithm, params.init=NULL, covar=NULL, diagnostic=FALSE, ...){
  ## Initialization for all algorithms
  algorithm <- match.arg(algorithm, choices=c("LEAP", "TS", "NTS", "THS"))
  fn <- function(x) {
    z <- -obj$fn(x)
    if(is.nan(z)){
      warning(paste("replacing NaN w/ Inf at:", paste(x, collapse=" ")))
      z <- Inf
    }
    return(z)
  }
  gr <- function(x) {
    z <- -as.vector(obj$gr(x))
    if(any(is.nan(z))){
      warning(paste("NaN gradient at:", paste(x, collapse=" ")))
      z <- rep(0, length(x))
    }
    return(z)
  }
  obj$env$beSilent()                  # silence console output
  ## argument checking
  if(is.null(params.init)){
    params.init <- obj$par
  } else if(length(params.init) != length(obj$par)){
    stop("params.init is wrong length")
  }
  ## Select and run the chain.
  if(algorithm == "LEAP")
    time <- system.time(mcmc.out <-
                          mcmc.leap(nsim=nsim, fn=fn, gr=gr, params.init=params.init, covar=covar,
                                    diagnostic=diagnostic, ...))
  else if(algorithm == "TS")
    time <- system.time(mcmc.out <-
                          mcmc.ts(nsim=nsim, fn=fn, gr=gr, params.init=params.init, covar=covar,
                                  diagnostic=diagnostic, ...))
  else if(algorithm == "NTS")
    time <- system.time(mcmc.out <-
                          mcmc.nts(nsim=nsim, fn=fn, gr=gr, params.init=params.init, covar=covar,
                                   diagnostic=diagnostic, ...))
  else if(algorithm == "THS")
    time <- system.time(mcmc.out <-
                          mcmc.ths(nsim=nsim, fn=fn, gr=gr, params.init=params.init, covar=covar,
                                   diagnostic=diagnostic, ...))
  ## Clean up returned output, a matrix if diag is FALSE, otherwise a list
  if(!diagnostic){
    mcmc.out <- as.data.frame(mcmc.out)
    names(mcmc.out) <- names(obj$par)
  } else {
    mcmc.out$time <- as.numeric(time[3])        # grab the elapsed time
    mcmc.out$par <- as.data.frame(mcmc.out$par)
    names(mcmc.out$par) <- names(obj$par)
  }
  return(invisible(mcmc.out))
}

#' [BETA VERSION] Draw MCMC samples from a model posterior using the
#' No-U-Turn (NUTS) sampler with dual averaging.
#'######## The LEAPFROG METHOD starts here
mcmc.leap <- function(nsim, fn, gr, params.init, max_doublings=4, eps=NULL, Madapt=NULL,
                      delta=0.5, covar=NULL, diagnostic=FALSE){
  ## If using covariance matrix and Cholesky decomposition, redefine
  ## these functions to include this transformation. The algorithm will
  ## work in the transformed space
  if(!is.null(covar)){
    fn2 <- function(theta) fn(chd %*% theta)
    gr2 <- function(theta) as.vector( t( gr(chd %*% theta) ) %*% chd )
    chd <- t(chol(covar))               # lower triangular Cholesky decomp.
    chd.inv <- solve(chd)               # inverse
    theta.cur <- chd.inv %*% params.init
  } else {
    fn2 <- fn; gr2 <- gr
    theta.cur <- params.init
  }
  theta.out <- matrix(NA, nrow=nsim, ncol=length(theta.cur))
  ## how many steps were taken at each iteration, useful for tuning
  j.results <- rep(NA, len=nsim)
  ## count the model calls; updated inside .buildtree. Some subtrees
  ## wont finish due to exit conditions so this is dynamic and not a
  ## simple formula like with HMC.
  info <- as.environment( list(n.calls = 0) )
  useDA <- is.null(eps)               # whether to use DA algorithm
  if(useDA){
    if(is.null(Madapt)){
      message("MCMC NUTS: Madapt not specified, defaulting to half of nsim")
      Madapt <- floor(nsim/2)
    }
    ## Initialize the dual-averaging algorithm.
    message(paste("MCMC NUTS: No eps given so using dual averaging during first", Madapt, "steps."))
    epsvec <- Hbar <- epsbar <- rep(NA, length=Madapt+1)
    eps <- epsvec[1] <- epsbar[1] <-
      .find.epsilon.leap(theta=theta.cur, fn=fn2, gr=gr2, eps=.1, verbose=FALSE)
    mu <- log(10*eps)
    Hbar[1] <- 0; gamma <- 0.05; t0 <- 10; kappa <- 0.75
  } else {
    ## dummy values to return
    epsvec <- epsbar <- Hbar <- NULL
  }
  ## Start of MCMC chain
  for(m in 1:nsim){
    ## initialize
    theta.out[m,] <- theta.minus <- theta.plus <- theta0 <- theta.cur
    r.cur <- r.plus <- r.minus <- r0 <- rnorm(length(theta.cur),0,1)
    ## Draw a slice variable u
    u <- .sample.u(theta=theta.cur, r=r.cur, fn=fn2)
    j <- 0; n <- 1; s <- 1
    while(s==1) {
      v <- sample(x=c(1,-1), size=1)
      if(v==1){
        ## move in right direction
        res <- .buildtree.leap(theta=theta.plus, r=r.plus, u=u, v=v,
                          j=j, eps=eps, theta0=theta0, r0=r0,
                          fn=fn2, gr=gr2, info=info)
        theta.plus <- res$theta.plus
        r.plus <- res$r.plus
      } else {
        ## move in left direction
        res <- .buildtree.leap(theta=theta.minus, r=r.minus, u=u, v=v,
                          j=j, eps=eps, theta0=theta0, r0=r0,
                          fn=fn2, gr=gr2, info=info)
        theta.minus <- res$theta.minus
        r.minus <- res$r.minus
      }
      ## test whether to accept this state
      if(is.na(res$s) | is.nan(res$s))  res$s <- 0
      if(res$s==1) {
        if(runif(n=1, min=0,max=1) <= res$n/n){
          theta.cur <- res$theta.prime
          theta.out[m,] <- res$theta.prime
        }
      }
      n <- n+res$n
      s <- res$s*.test.nuts(theta.plus, theta.minus, r.plus, r.minus)
      ## Stop trajectory if there are any problems, probably happens
      ## when jumping way too far into the tails and the model isn't
      ## defined
      if(is.na(s) | is.nan(s))  s <- 0
      j <- j+1
      ## Stop doubling if too many or it's diverged enough
      if(j>max_doublings & s) {
        ## warning("j larger than max_doublings, skipping to next m")
        break
      }
    }
    j.results[m] <- j-1
    if(useDA){
      ## Do the adapting of eps.
      if(m <= Madapt){
        Hbar[m+1] <- (1-1/(m+t0))*Hbar[m] +
          (delta-res$alpha/res$nalpha)/(m+t0)
        ## If logalpha not defined, skip this updating step and use
        ## the last one.
        if(is.nan(Hbar[m+1])) Hbar[m+1] <- abs(Hbar[m])
        
        logeps <- mu-sqrt(m)*Hbar[m+1]/gamma
        epsvec[m+1] <- exp(logeps)
        logepsbar <- m^(-kappa)*logeps + (1-m^(-kappa))*log(epsbar[m])
        epsbar[m+1] <- exp(logepsbar)
        eps <- epsvec[m+1]
      } else {
        eps <- epsbar[Madapt]*runif(1,.9,1.1)
      }
    }
  } ## end of MCMC loop
  ## Back transform parameters if covar is used
  if(!is.null(covar)) {
    theta.out <- t(apply(theta.out, 1, function(x) chd %*% x))
  }
  j.stats <- 2^(c(min(j.results), median(j.results), max(j.results)))
  if(useDA)
    message(paste("MCMC NUTS: Dual averaging final average eps =", round(epsbar[Madapt], 3)))
  message(paste0("MCMC NUTS: Approximate leapfrog steps(min, median, max)=(",
                 paste(j.stats, collapse=","), ")"))
  if(diagnostic){
    return(list(par=theta.out, steps.taken= 2^j.results,
                n.calls=info$n.calls, epsvec=epsvec, epsbar=epsbar, Hbar=Hbar))
  } else {
    return(theta.out)
  }
}

# Draw a slice sample for given position and momentum variables
.sample.u <- function(theta, r, fn)
  runif(n=1, min=0, max=exp(.calculate.H(theta=theta,r=r, fn=fn)))
# Calculate the Hamiltonian value for position and momentum variables.
#
# @details This function currently assumes iid standard normal momentum
# variables.
.calculate.H <- function(theta, r, fn) fn(theta)-(1/2)*sum(r^2)
# Test whether a "U-turn" has occured in a branch of the binary tree
# created by \ref\code{.buildtree.eps} function.
.test.nuts <- function(theta.plus, theta.minus, r.plus, r.minus){
  theta.temp <- (theta.plus-theta.minus)
  as.numeric( theta.temp %*% r.minus >= 0 | theta.temp %*% r.plus >= 0)
}

# A recursive function that builds a leapfrog trajectory using a balanced
# binary tree.
#
# @references This is from the No-U-Turn sampler with dual averaging
# (algorithm 6) of Hoffman and Gelman (2014).
#
# @details The function repeatedly doubles (in a random direction) until
# either a U-turn occurs or the trajectory becomes unstable. This is the
# 'efficient' version that samples uniformly from the path without storing
# it. Thus the function returns a single proposed value and not the whole
# trajectory.
.buildtree.leap <- function(theta, r, u, v, j, eps, theta0, r0, fn, gr,
                       delta.max=1000, info = environment() ){
  if(j==0){
    ## base case, take one step in direction v ##### THIS IS WHERE IT IS LEAPFROG SPECIFIC ####
    eps <- v*eps
    r <- r+(eps/2)*gr(theta)
    theta <- theta+eps*r
    r <- r+(eps/2)*gr(theta)
    ## verify valid trajectory
    H <- .calculate.H(theta=theta, r=r, fn=fn)
    s <- H-log(u) + delta.max > 0
    if(is.na(s) | is.nan(s)) s <- 0
    n <- log(u) <= H
    ## ## Useful code for debugging. Returns entire path to global env.
    ## if(!exists('theta.trajectory'))
    ##     theta.trajectory <<- theta
    ## else
    ##     theta.trajectory <<- rbind(theta.trajectory, theta)
    temp <- .calculate.H(theta=theta, r=r, fn=fn)-
      .calculate.H(theta=theta0, r=r0, fn=fn)
    alpha <- min(exp(temp),1)
    info$n.calls <- info$n.calls + 5
    return(list(theta.minus=theta, theta.plus=theta, theta.prime=theta, r.minus=r,
                r.plus=r, s=s, n=n, alpha=alpha, nalpha=1))
  } else {
    ## recursion - build left and right subtrees
    xx <- .buildtree.leap(theta=theta, r=r, u=u, v=v, j=j-1, eps=eps,
                     theta0=theta0, r0=r0, fn=fn, gr=gr, info=info)
    theta.minus <- xx$theta.minus
    theta.plus <- xx$theta.plus
    theta.prime <- xx$theta.prime
    r.minus <- xx$r.minus
    r.plus <- xx$r.plus
    alpha <- xx$alpha
    nalpha <- xx$nalpha
    s <- xx$s
    if(is.na(s) | is.nan(s)) s <- 0
    nprime <- xx$n
    ## If it didn't fail, update the above quantities
    if(s==1){
      if(v== -1){
        yy <- .buildtree.leap(theta=theta.minus, r=r.minus, u=u, v=v,
                         j=j-1, eps=eps, theta0=theta0, r0=r0,
                         fn=fn, gr=gr, info=info)
        theta.minus <- yy$theta.minus
        r.minus <- yy$r.minus
      } else {
        yy <- .buildtree.leap(theta=theta.plus, r=r.plus, u=u, v=v,
                         j=j-1, eps=eps, theta0=theta0, r0=r0,
                         fn=fn, gr=gr, info=info)
        theta.plus <- yy$theta.plus
        r.plus <- yy$r.plus
      }
      ## This isn't in the paper but if both slice variables failed,
      ## then you get 0/0. So I skip this test. Likewise if model
      ## throwing errors, don't keep that theta.
      nprime <- yy$n+ xx$n
      if(!is.finite(nprime)) nprime <- 0
      if(nprime!=0){
        ## choose whether to keep this theta
        if(runif(n=1, min=0, max=1) <= yy$n/nprime)
          theta.prime <- yy$theta.prime
      }
      ## check for valid proposal
      test <- .test.nuts(theta.plus=theta.plus,
                         theta.minus=theta.minus, r.plus=r.plus,
                         r.minus=r.minus)
      ## if(!test) warning(paste("U turn at j=", j))
      ## check if any of the stopping conditions were met
      s <- xx$s*yy$s*test
    }
    return(list(theta.minus=theta.minus, theta.plus=theta.plus,
                theta.prime=theta.prime,
                r.minus=r.minus, r.plus=r.plus, s=s, n=nprime,
                alpha=alpha, nalpha=1))
  }
}

# Estimate a reasonable starting value for epsilon (step size) for a given
# model, for use with Hamiltonian MCMC algorithms.
#
# This is Algorithm 4 from Hoffman and Gelman (2010) and is used in the
# dual-averaging algorithms for both HMC and NUTS to find a reasonable
# starting value.
# @title Estimate step size for Hamiltonian MCMC algorithms
# @param theta An initial parameter vector.
# @param fn A function returning the log-likelihood (not the negative of
# it) for a given parameter vector.
# @param gr A function returning the gradient of the log-likelihood of a
# model.
# @param eps A value for espilon to initiate the algorithm. Defaults to
# 1. If this is far too big the algorithm won't work well and an
# alternative value can be used.
# @return Returns the "reasonable" espilon invisible, while printing how
# many steps to reach it.
# @details The algorithm uses a while loop and will break after 50
# iterations.
#
.find.epsilon.leap <- function(theta,  fn, gr, eps=1, verbose=TRUE){
  r <- rnorm(n=length(theta), mean=0, sd=1)
  ## Do one leapfrog step
  r.new <- r+(eps/2)*gr(theta)
  theta.new <- theta+eps*r.new
  r.new <- r.new+(eps/2)*gr(theta.new)
  H1 <- .calculate.H(theta=theta, r=r, fn=fn)
  H2 <- .calculate.H(theta=theta.new, r=r.new, fn=fn)
  a <- 2*(exp(H2)/exp(H1)>.5)-1
  ## If jumped into bad region, a can be NaN so setup algorithm to keep
  ## halving eps instead of throwing error
  if(!is.finite(a)) a <- -1
  k <- 1
  ## Similarly, keep going if there are infinite values
  while (!is.finite(H1) | !is.finite(H2) | a*H2-a*H1 > -a*log(2)) {
    eps <- (2^a)*eps
    ## Do one leapfrog step
    r.new <- r+(eps/2)*gr(theta)
    theta.new <- theta+eps*r.new
    r.new <- r.new+(eps/2)*gr(theta.new)
    H2 <- .calculate.H(theta=theta.new, r=r.new, fn=fn)
    k <- k+1
    if(k>50) {
      stop("More than 50 iterations to find reasonable eps. Model is likely misspecified or some other issue.")
    }
  }
  if(verbose) message(paste("Reasonable epsilon=", eps, "found after", k, "steps"))
  return(invisible(eps))
}

## Copyright (C) 2015 Cole Monnahan
## License: GPL-2

#' [BETA VERSION] Draw MCMC samples from a model posterior using the
#' No-U-Turn (NUTS) sampler with dual averaging and TWO STAGE METHOD
mcmc.ts <- function(nsim, fn, gr, params.init, max_doublings=4, eps=NULL, Madapt=NULL,
                    delta=0.5, covar=NULL, diagnostic=FALSE){
  ## If using covariance matrix and Cholesky decomposition, redefine
  ## these functions to include this transformation. The algorithm will
  ## work in the transformed space
  if(!is.null(covar)){
    fn2 <- function(theta) fn(chd %*% theta)
    gr2 <- function(theta) as.vector( t( gr(chd %*% theta) ) %*% chd )
    chd <- t(chol(covar))               # lower triangular Cholesky decomp.
    chd.inv <- solve(chd)               # inverse
    theta.cur <- chd.inv %*% params.init
  } else {
    fn2 <- fn; gr2 <- gr
    theta.cur <- params.init
  }
  theta.out <- matrix(NA, nrow=nsim, ncol=length(theta.cur))
  ## how many steps were taken at each iteration, useful for tuning
  j.results <- rep(NA, len=nsim)
  ## count the model calls; updated inside .buildtree. Some subtrees
  ## wont finish due to exit conditions so this is dynamic and not a
  ## simple formula like with HMC.
  info <- as.environment( list(n.calls = 0) )
  useDA <- is.null(eps)               # whether to use DA algorithm
  if(useDA){
    if(is.null(Madapt)){
      message("MCMC NUTS: Madapt not specified, defaulting to half of nsim")
      Madapt <- floor(nsim/2)
    }
    ## Initialize the dual-averaging algorithm.
    message(paste("MCMC NUTS: No eps given so using dual averaging during first", Madapt, "steps."))
    epsvec <- Hbar <- epsbar <- rep(NA, length=Madapt+1)
    eps <- epsvec[1] <- epsbar[1] <-
      .find.epsilon.ts(theta=theta.cur, fn=fn2, gr=gr2, eps=.1, verbose=FALSE)
    mu <- log(10*eps)
    Hbar[1] <- 0; gamma <- 0.05; t0 <- 10; kappa <- 0.75
  } else {
    ## dummy values to return
    epsvec <- epsbar <- Hbar <- NULL
  }
  ## Start of MCMC chain
  for(m in 1:nsim){
    ## initialize
    theta.out[m,] <- theta.minus <- theta.plus <- theta0 <- theta.cur
    r.cur <- r.plus <- r.minus <- r0 <- rnorm(length(theta.cur),0,1)
    ## Draw a slice variable u
    u <- .sample.u(theta=theta.cur, r=r.cur, fn=fn2)
    j <- 0; n <- 1; s <- 1
    while(s==1) {
      v <- sample(x=c(1,-1), size=1)
      if(v==1){
        ## move in right direction
        res <- .buildtree.ts(theta=theta.plus, r=r.plus, u=u, v=v,
                          j=j, eps=eps, theta0=theta0, r0=r0,
                          fn=fn2, gr=gr2, info=info)
        theta.plus <- res$theta.plus
        r.plus <- res$r.plus
      } else {
        ## move in left direction
        res <- .buildtree.ts(theta=theta.minus, r=r.minus, u=u, v=v,
                          j=j, eps=eps, theta0=theta0, r0=r0,
                          fn=fn2, gr=gr2, info=info)
        theta.minus <- res$theta.minus
        r.minus <- res$r.minus
      }
      ## test whether to accept this state
      if(is.na(res$s) | is.nan(res$s))  res$s <- 0
      if(res$s==1) {
        if(runif(n=1, min=0,max=1) <= res$n/n){
          theta.cur <- res$theta.prime
          theta.out[m,] <- res$theta.prime
        }
      }
      n <- n+res$n
      s <- res$s*.test.nuts(theta.plus, theta.minus, r.plus, r.minus)
      ## Stop trajectory if there are any problems, probably happens
      ## when jumping way too far into the tails and the model isn't
      ## defined
      if(is.na(s) | is.nan(s))  s <- 0
      j <- j+1
      ## Stop doubling if too many or it's diverged enough
      if(j>max_doublings & s) {
        ## warning("j larger than max_doublings, skipping to next m")
        break
      }
    }
    j.results[m] <- j-1
    if(useDA){
      ## Do the adapting of eps.
      if(m <= Madapt){
        Hbar[m+1] <- (1-1/(m+t0))*Hbar[m] +
          (delta-res$alpha/res$nalpha)/(m+t0)
        ## If logalpha not defined, skip this updating step and use
        ## the last one.
        if(is.nan(Hbar[m+1])) Hbar[m+1] <- abs(Hbar[m])
        
        logeps <- mu-sqrt(m)*Hbar[m+1]/gamma
        epsvec[m+1] <- exp(logeps)
        logepsbar <- m^(-kappa)*logeps + (1-m^(-kappa))*log(epsbar[m])
        epsbar[m+1] <- exp(logepsbar)
        eps <- epsvec[m+1]
      } else {
        eps <- epsbar[Madapt]*runif(1,.9,1.1)
      }
    }
  } ## end of MCMC loop
  ## Back transform parameters if covar is used
  if(!is.null(covar)) {
    theta.out <- t(apply(theta.out, 1, function(x) chd %*% x))
  }
  j.stats <- 2^(c(min(j.results), median(j.results), max(j.results)))
  if(useDA)
    message(paste("MCMC NUTS: Dual averaging final average eps =", round(epsbar[Madapt], 3)))
  message(paste0("MCMC NUTS: Approximate two stage steps(min, median, max)=(",
                 paste(j.stats, collapse=","), ")"))
  if(diagnostic){
    return(list(par=theta.out, steps.taken= 2^j.results,
                n.calls=info$n.calls, epsvec=epsvec, epsbar=epsbar, Hbar=Hbar))
  } else {
    return(theta.out)
  }
}


#### THIS IS WHERE TWO STAGE SPEFICIC
.buildtree.ts <- function(theta, r, u, v, j, eps, theta0, r0, fn, gr,
                       delta.max=1000, info = environment() ){
  if(j==0){
    ## base case, take one step in direction v
    eps <- v*eps 
    #### two stage method ####
    a1 <- 0.21132
    theta <- theta + a1*eps*r
    r <- r + (eps/2)*gr(theta)
    theta <- theta + (1-2*a1)*eps*r
    r <- r + (eps/2)* gr(theta)
    theta <- theta + a1*eps*r
    
    ## verify valid trajectory
    H <- .calculate.H(theta=theta, r=r, fn=fn)
    s <- H-log(u) + delta.max > 0
    if(is.na(s) | is.nan(s)) s <- 0
    n <- log(u) <= H
    ## ## Useful code for debugging. Returns entire path to global env.
    ## if(!exists('theta.trajectory'))
    ##     theta.trajectory <<- theta
    ## else
    ##     theta.trajectory <<- rbind(theta.trajectory, theta)
    temp <- .calculate.H(theta=theta, r=r, fn=fn)-
      .calculate.H(theta=theta0, r=r0, fn=fn)
    alpha <- min(exp(temp),1)
    info$n.calls <- info$n.calls + 5
    return(list(theta.minus=theta, theta.plus=theta, theta.prime=theta, r.minus=r,
                r.plus=r, s=s, n=n, alpha=alpha, nalpha=1))
  } else {
    ## recursion - build left and right subtrees
    xx <- .buildtree.ts(theta=theta, r=r, u=u, v=v, j=j-1, eps=eps,
                     theta0=theta0, r0=r0, fn=fn, gr=gr, info=info)
    theta.minus <- xx$theta.minus
    theta.plus <- xx$theta.plus
    theta.prime <- xx$theta.prime
    r.minus <- xx$r.minus
    r.plus <- xx$r.plus
    alpha <- xx$alpha
    nalpha <- xx$nalpha
    s <- xx$s
    if(is.na(s) | is.nan(s)) s <- 0
    nprime <- xx$n
    ## If it didn't fail, update the above quantities
    if(s==1){
      if(v== -1){
        yy <- .buildtree.ts(theta=theta.minus, r=r.minus, u=u, v=v,
                         j=j-1, eps=eps, theta0=theta0, r0=r0,
                         fn=fn, gr=gr, info=info)
        theta.minus <- yy$theta.minus
        r.minus <- yy$r.minus
      } else {
        yy <- .buildtree.ts(theta=theta.plus, r=r.plus, u=u, v=v,
                         j=j-1, eps=eps, theta0=theta0, r0=r0,
                         fn=fn, gr=gr, info=info)
        theta.plus <- yy$theta.plus
        r.plus <- yy$r.plus
      }
      ## This isn't in the paper but if both slice variables failed,
      ## then you get 0/0. So I skip this test. Likewise if model
      ## throwing errors, don't keep that theta.
      nprime <- yy$n+ xx$n
      if(!is.finite(nprime)) nprime <- 0
      if(nprime!=0){
        ## choose whether to keep this theta
        if(runif(n=1, min=0, max=1) <= yy$n/nprime)
          theta.prime <- yy$theta.prime
      }
      ## check for valid proposal
      test <- .test.nuts(theta.plus=theta.plus,
                         theta.minus=theta.minus, r.plus=r.plus,
                         r.minus=r.minus)
      ## if(!test) warning(paste("U turn at j=", j))
      ## check if any of the stopping conditions were met
      s <- xx$s*yy$s*test
    }
    return(list(theta.minus=theta.minus, theta.plus=theta.plus,
                theta.prime=theta.prime,
                r.minus=r.minus, r.plus=r.plus, s=s, n=nprime,
                alpha=alpha, nalpha=1))
  }
}

#Finding epsilon
.find.epsilon.ts <- function(theta,  fn, gr, eps=1, verbose=TRUE){
  r <- rnorm(n=length(theta), mean=0, sd=1)
  ## Do one two stage step
  a1 <- 0.21132
  theta.new <- theta + a1*eps*r
  r.new <- r + (eps/2)*gr(theta.new)
  theta.new <- theta.new + (1-2*a1)*eps*r.new
  r.new <- r.new + (eps/2)* gr(theta.new)
  theta.new <- theta.new + a1*eps*r.new
  
  H1 <- .calculate.H(theta=theta, r=r, fn=fn)
  H2 <- .calculate.H(theta=theta.new, r=r.new, fn=fn)
  a <- 2*(exp(H2)/exp(H1)>.5)-1
  ## If jumped into bad region, a can be NaN so setup algorithm to keep
  ## halving eps instead of throwing error
  if(!is.finite(a)) a <- -1
  k <- 1
  ## Similarly, keep going if there are infinite values
  while (!is.finite(H1) | !is.finite(H2) | a*H2-a*H1 > -a*log(2)) {
    eps <- (2^a)*eps
    ## Do one leapfrog step
    a1 <- 0.21132
    theta.new <- theta + a1*eps*r
    r.new <- r + (eps/2)*gr(theta.new)
    theta.new <- theta.new + (1-2*a1)*eps*r.new
    r.new <- r.new + (eps/2)* gr(theta.new)
    theta.new <- theta.new + a1*eps*r.new
    H2 <- .calculate.H(theta=theta.new, r=r.new, fn=fn)
    k <- k+1
    if(k>50) {
      stop("More than 50 iterations to find reasonable eps. Model is likely misspecified or some other issue.")
    }
  }
  if(verbose) message(paste("Reasonable epsilon=", eps, "found after", k, "steps"))
  return(invisible(eps))
}


## Copyright (C) 2015 Cole Monnahan
## License: GPL-2

#' [BETA VERSION] Draw MCMC samples from a model posterior using the
#' No-U-Turn (NUTS) sampler with dual averaging and NEW TWO STAGE METHOD
mcmc.nts <- function(nsim, fn, gr, params.init, max_doublings=4, eps=NULL, Madapt=NULL,
                     delta=0.5, covar=NULL, diagnostic=FALSE){
  ## If using covariance matrix and Cholesky decomposition, redefine
  ## these functions to include this transformation. The algorithm will
  ## work in the transformed space
  if(!is.null(covar)){
    fn2 <- function(theta) fn(chd %*% theta)
    gr2 <- function(theta) as.vector( t( gr(chd %*% theta) ) %*% chd )
    chd <- t(chol(covar))               # lower triangular Cholesky decomp.
    chd.inv <- solve(chd)               # inverse
    theta.cur <- chd.inv %*% params.init
  } else {
    fn2 <- fn; gr2 <- gr
    theta.cur <- params.init
  }
  theta.out <- matrix(NA, nrow=nsim, ncol=length(theta.cur))
  ## how many steps were taken at each iteration, useful for tuning
  j.results <- rep(NA, len=nsim)
  ## count the model calls; updated inside .buildtree. Some subtrees
  ## wont finish due to exit conditions so this is dynamic and not a
  ## simple formula like with HMC.
  info <- as.environment( list(n.calls = 0) )
  useDA <- is.null(eps)               # whether to use DA algorithm
  if(useDA){
    if(is.null(Madapt)){
      message("MCMC NUTS: Madapt not specified, defaulting to half of nsim")
      Madapt <- floor(nsim/2)
    }
    ## Initialize the dual-averaging algorithm.
    message(paste("MCMC NUTS: No eps given so using dual averaging during first", Madapt, "steps."))
    epsvec <- Hbar <- epsbar <- rep(NA, length=Madapt+1)
    eps <- epsvec[1] <- epsbar[1] <-
      .find.epsilon.nts(theta=theta.cur, fn=fn2, gr=gr2, eps=.1, verbose=FALSE)
    mu <- log(10*eps)
    Hbar[1] <- 0; gamma <- 0.05; t0 <- 10; kappa <- 0.75
  } else {
    ## dummy values to return
    epsvec <- epsbar <- Hbar <- NULL
  }
  ## Start of MCMC chain
  for(m in 1:nsim){
    ## initialize
    theta.out[m,] <- theta.minus <- theta.plus <- theta0 <- theta.cur
    r.cur <- r.plus <- r.minus <- r0 <- rnorm(length(theta.cur),0,1)
    ## Draw a slice variable u
    u <- .sample.u(theta=theta.cur, r=r.cur, fn=fn2)
    j <- 0; n <- 1; s <- 1
    while(s==1) {
      v <- sample(x=c(1,-1), size=1)
      if(v==1){
        ## move in right direction
        res <- .buildtree.nts(theta=theta.plus, r=r.plus, u=u, v=v,
                          j=j, eps=eps, theta0=theta0, r0=r0,
                          fn=fn2, gr=gr2, info=info)
        theta.plus <- res$theta.plus
        r.plus <- res$r.plus
      } else {
        ## move in left direction
        res <- .buildtree.nts(theta=theta.minus, r=r.minus, u=u, v=v,
                          j=j, eps=eps, theta0=theta0, r0=r0,
                          fn=fn2, gr=gr2, info=info)
        theta.minus <- res$theta.minus
        r.minus <- res$r.minus
      }
      ## test whether to accept this state
      if(is.na(res$s) | is.nan(res$s))  res$s <- 0
      if(res$s==1) {
        if(runif(n=1, min=0,max=1) <= res$n/n){
          theta.cur <- res$theta.prime
          theta.out[m,] <- res$theta.prime
        }
      }
      n <- n+res$n
      s <- res$s*.test.nuts(theta.plus, theta.minus, r.plus, r.minus)
      ## Stop trajectory if there are any problems, probably happens
      ## when jumping way too far into the tails and the model isn't
      ## defined
      if(is.na(s) | is.nan(s))  s <- 0
      j <- j+1
      ## Stop doubling if too many or it's diverged enough
      if(j>max_doublings & s) {
        ## warning("j larger than max_doublings, skipping to next m")
        break
      }
    }
    j.results[m] <- j-1
    if(useDA){
      ## Do the adapting of eps.
      if(m <= Madapt){
        Hbar[m+1] <- (1-1/(m+t0))*Hbar[m] +
          (delta-res$alpha/res$nalpha)/(m+t0)
        ## If logalpha not defined, skip this updating step and use
        ## the last one.
        if(is.nan(Hbar[m+1])) Hbar[m+1] <- abs(Hbar[m])
        
        logeps <- mu-sqrt(m)*Hbar[m+1]/gamma
        epsvec[m+1] <- exp(logeps)
        logepsbar <- m^(-kappa)*logeps + (1-m^(-kappa))*log(epsbar[m])
        epsbar[m+1] <- exp(logepsbar)
        eps <- epsvec[m+1]
      } else {
        eps <- epsbar[Madapt]*runif(1,.9,1.1)
      }
    }
  } ## end of MCMC loop
  ## Back transform parameters if covar is used
  if(!is.null(covar)) {
    theta.out <- t(apply(theta.out, 1, function(x) chd %*% x))
  }
  j.stats <- 2^(c(min(j.results), median(j.results), max(j.results)))
  if(useDA)
    message(paste("MCMC NUTS: Dual averaging final average eps =", round(epsbar[Madapt], 3)))
  message(paste0("MCMC NUTS: Approximate new two stage steps(min, median, max)=(",
                 paste(j.stats, collapse=","), ")"))
  if(diagnostic){
    return(list(par=theta.out, steps.taken= 2^j.results,
                n.calls=info$n.calls, epsvec=epsvec, epsbar=epsbar, Hbar=Hbar))
  } else {
    return(theta.out)
  }
}

#### THis is NEW TWO STAGE SPECIFIC ####
.buildtree.nts <- function(theta, r, u, v, j, eps, theta0, r0, fn, gr,
                       delta.max=1000, info = environment() ){
  if(j==0){
    ## base case, take one step in direction v
    eps <- v*eps 
    #### New two stage method ####
    a1 <- 0.19098
    theta <- theta + a1*eps*r
    r <- r + (eps/2)*gr(theta)
    theta <- theta + (1-2*a1)*eps*r
    r <- r + (eps/2)* gr(theta)
    theta <- theta + a1*eps*r
    
    ## verify valid trajectory
    H <- .calculate.H(theta=theta, r=r, fn=fn)
    s <- H-log(u) + delta.max > 0
    if(is.na(s) | is.nan(s)) s <- 0
    n <- log(u) <= H
    ## ## Useful code for debugging. Returns entire path to global env.
    ## if(!exists('theta.trajectory'))
    ##     theta.trajectory <<- theta
    ## else
    ##     theta.trajectory <<- rbind(theta.trajectory, theta)
    temp <- .calculate.H(theta=theta, r=r, fn=fn)-
      .calculate.H(theta=theta0, r=r0, fn=fn)
    alpha <- min(exp(temp),1)
    info$n.calls <- info$n.calls + 5
    return(list(theta.minus=theta, theta.plus=theta, theta.prime=theta, r.minus=r,
                r.plus=r, s=s, n=n, alpha=alpha, nalpha=1))
  } else {
    ## recursion - build left and right subtrees
    xx <- .buildtree.nts(theta=theta, r=r, u=u, v=v, j=j-1, eps=eps,
                     theta0=theta0, r0=r0, fn=fn, gr=gr, info=info)
    theta.minus <- xx$theta.minus
    theta.plus <- xx$theta.plus
    theta.prime <- xx$theta.prime
    r.minus <- xx$r.minus
    r.plus <- xx$r.plus
    alpha <- xx$alpha
    nalpha <- xx$nalpha
    s <- xx$s
    if(is.na(s) | is.nan(s)) s <- 0
    nprime <- xx$n
    ## If it didn't fail, update the above quantities
    if(s==1){
      if(v== -1){
        yy <- .buildtree.nts(theta=theta.minus, r=r.minus, u=u, v=v,
                         j=j-1, eps=eps, theta0=theta0, r0=r0,
                         fn=fn, gr=gr, info=info)
        theta.minus <- yy$theta.minus
        r.minus <- yy$r.minus
      } else {
        yy <- .buildtree.nts(theta=theta.plus, r=r.plus, u=u, v=v,
                         j=j-1, eps=eps, theta0=theta0, r0=r0,
                         fn=fn, gr=gr, info=info)
        theta.plus <- yy$theta.plus
        r.plus <- yy$r.plus
      }
      ## This isn't in the paper but if both slice variables failed,
      ## then you get 0/0. So I skip this test. Likewise if model
      ## throwing errors, don't keep that theta.
      nprime <- yy$n+ xx$n
      if(!is.finite(nprime)) nprime <- 0
      if(nprime!=0){
        ## choose whether to keep this theta
        if(runif(n=1, min=0, max=1) <= yy$n/nprime)
          theta.prime <- yy$theta.prime
      }
      ## check for valid proposal
      test <- .test.nuts(theta.plus=theta.plus,
                         theta.minus=theta.minus, r.plus=r.plus,
                         r.minus=r.minus)
      ## if(!test) warning(paste("U turn at j=", j))
      ## check if any of the stopping conditions were met
      s <- xx$s*yy$s*test
    }
    return(list(theta.minus=theta.minus, theta.plus=theta.plus,
                theta.prime=theta.prime,
                r.minus=r.minus, r.plus=r.plus, s=s, n=nprime,
                alpha=alpha, nalpha=1))
  }
}

#Find epsilon
.find.epsilon.nts <- function(theta,  fn, gr, eps=1, verbose=TRUE){
  r <- rnorm(n=length(theta), mean=0, sd=1)
  ## Do one NEW TWO stage step
  a1 <- 0.19098
  theta.new <- theta + a1*eps*r
  r.new <- r + (eps/2)*gr(theta.new)
  theta.new <- theta.new + (1-2*a1)*eps*r.new
  r.new <- r.new + (eps/2)* gr(theta.new)
  theta.new <- theta.new + a1*eps*r.new
  
  H1 <- .calculate.H(theta=theta, r=r, fn=fn)
  H2 <- .calculate.H(theta=theta.new, r=r.new, fn=fn)
  a <- 2*(exp(H2)/exp(H1)>.5)-1
  ## If jumped into bad region, a can be NaN so setup algorithm to keep
  ## halving eps instead of throwing error
  if(!is.finite(a)) a <- -1
  k <- 1
  ## Similarly, keep going if there are infinite values
  while (!is.finite(H1) | !is.finite(H2) | a*H2-a*H1 > -a*log(2)) {
    eps <- (2^a)*eps
    ## Do one leapfrog step
    a1 <- 0.19098
    theta.new <- theta + a1*eps*r
    r.new <- r + (eps/2)*gr(theta.new)
    theta.new <- theta.new + (1-2*a1)*eps*r.new
    r.new <- r.new + (eps/2)* gr(theta.new)
    theta.new <- theta.new + a1*eps*r.new
    H2 <- .calculate.H(theta=theta.new, r=r.new, fn=fn)
    k <- k+1
    if(k>50) {
      stop("More than 50 iterations to find reasonable eps. Model is likely misspecified or some other issue.")
    }
  }
  if(verbose) message(paste("Reasonable epsilon=", eps, "found after", k, "steps"))
  return(invisible(eps))
}


## Copyright (C) 2015 Cole Monnahan
## License: GPL-2

#' [BETA VERSION] Draw MCMC samples from a model posterior using the
#' No-U-Turn (NUTS) sampler with dual averaging and THREE STAGE METHOD.
mcmc.ths <- function(nsim, fn, gr, params.init, max_doublings=4, eps=NULL, Madapt=NULL,
                     delta=0.5, covar=NULL, diagnostic=FALSE){
  ## If using covariance matrix and Cholesky decomposition, redefine
  ## these functions to include this transformation. The algorithm will
  ## work in the transformed space
  if(!is.null(covar)){
    fn2 <- function(theta) fn(chd %*% theta)
    gr2 <- function(theta) as.vector( t( gr(chd %*% theta) ) %*% chd )
    chd <- t(chol(covar))               # lower triangular Cholesky decomp.
    chd.inv <- solve(chd)               # inverse
    theta.cur <- chd.inv %*% params.init
  } else {
    fn2 <- fn; gr2 <- gr
    theta.cur <- params.init
  }
  theta.out <- matrix(NA, nrow=nsim, ncol=length(theta.cur))
  ## how many steps were taken at each iteration, useful for tuning
  j.results <- rep(NA, len=nsim)
  ## count the model calls; updated inside .buildtree. Some subtrees
  ## wont finish due to exit conditions so this is dynamic and not a
  ## simple formula like with HMC.
  info <- as.environment( list(n.calls = 0) )
  useDA <- is.null(eps)               # whether to use DA algorithm
  if(useDA){
    if(is.null(Madapt)){
      message("MCMC NUTS: Madapt not specified, defaulting to half of nsim")
      Madapt <- floor(nsim/2)
    }
    ## Initialize the dual-averaging algorithm.
    message(paste("MCMC NUTS: No eps given so using dual averaging during first", Madapt, "steps."))
    epsvec <- Hbar <- epsbar <- rep(NA, length=Madapt+1)
    eps <- epsvec[1] <- epsbar[1] <-
      .find.epsilon.ths(theta=theta.cur, fn=fn2, gr=gr2, eps=.1, verbose=FALSE)
    mu <- log(10*eps)
    Hbar[1] <- 0; gamma <- 0.05; t0 <- 10; kappa <- 0.75
  } else {
    ## dummy values to return
    epsvec <- epsbar <- Hbar <- NULL
  }
  ## Start of MCMC chain
  for(m in 1:nsim){
    ## initialize
    theta.out[m,] <- theta.minus <- theta.plus <- theta0 <- theta.cur
    r.cur <- r.plus <- r.minus <- r0 <- rnorm(length(theta.cur),0,1)
    ## Draw a slice variable u
    u <- .sample.u(theta=theta.cur, r=r.cur, fn=fn2)
    j <- 0; n <- 1; s <- 1
    while(s==1) {
      v <- sample(x=c(1,-1), size=1)
      if(v==1){
        ## move in right direction
        res <- .buildtree.ths(theta=theta.plus, r=r.plus, u=u, v=v,
                          j=j, eps=eps, theta0=theta0, r0=r0,
                          fn=fn2, gr=gr2, info=info)
        theta.plus <- res$theta.plus
        r.plus <- res$r.plus
      } else {
        ## move in left direction
        res <- .buildtree.ths(theta=theta.minus, r=r.minus, u=u, v=v,
                          j=j, eps=eps, theta0=theta0, r0=r0,
                          fn=fn2, gr=gr2, info=info)
        theta.minus <- res$theta.minus
        r.minus <- res$r.minus
      }
      ## test whether to accept this state
      if(is.na(res$s) | is.nan(res$s))  res$s <- 0
      if(res$s==1) {
        if(runif(n=1, min=0,max=1) <= res$n/n){
          theta.cur <- res$theta.prime
          theta.out[m,] <- res$theta.prime
        }
      }
      n <- n+res$n
      s <- res$s*.test.nuts(theta.plus, theta.minus, r.plus, r.minus)
      ## Stop trajectory if there are any problems, probably happens
      ## when jumping way too far into the tails and the model isn't
      ## defined
      if(is.na(s) | is.nan(s))  s <- 0
      j <- j+1
      ## Stop doubling if too many or it's diverged enough
      if(j>max_doublings & s) {
        ## warning("j larger than max_doublings, skipping to next m")
        break
      }
    }
    j.results[m] <- j-1
    if(useDA){
      ## Do the adapting of eps.
      if(m <= Madapt){
        Hbar[m+1] <- (1-1/(m+t0))*Hbar[m] +
          (delta-res$alpha/res$nalpha)/(m+t0)
        ## If logalpha not defined, skip this updating step and use
        ## the last one.
        if(is.nan(Hbar[m+1])) Hbar[m+1] <- abs(Hbar[m])
        
        logeps <- mu-sqrt(m)*Hbar[m+1]/gamma
        epsvec[m+1] <- exp(logeps)
        logepsbar <- m^(-kappa)*logeps + (1-m^(-kappa))*log(epsbar[m])
        epsbar[m+1] <- exp(logepsbar)
        eps <- epsvec[m+1]
      } else {
        eps <- epsbar[Madapt]*runif(1,.9,1.1)
      }
    }
  } ## end of MCMC loop
  ## Back transform parameters if covar is used
  if(!is.null(covar)) {
    theta.out <- t(apply(theta.out, 1, function(x) chd %*% x))
  }
  j.stats <- 2^(c(min(j.results), median(j.results), max(j.results)))
  if(useDA)
    message(paste("MCMC NUTS: Dual averaging final average eps =", round(epsbar[Madapt], 3)))
  message(paste0("MCMC NUTS: Approximate three stage steps(min, median, max)=(",
                 paste(j.stats, collapse=","), ")"))
  if(diagnostic){
    return(list(par=theta.out, steps.taken= 2^j.results,
                n.calls=info$n.calls, epsvec=epsvec, epsbar=epsbar, Hbar=Hbar))
  } else {
    return(theta.out)
  }
}

##### THIS IS THREE STAGE SPECIFIC ####
.buildtree.ths <- function(theta, r, u, v, j, eps, theta0, r0, fn, gr,
                       delta.max=1000, info = environment() ){
  if(j==0){
    ## base case, take one step in direction v
    eps <- v*eps
    #### three stage ####
    a1 <- 0.11888
    b1 <- 0.29620
    theta <- theta + a1*eps*r
    r <- r + b1*eps*gr(theta)
    theta <- theta + (0.5 - a1)*eps*r
    r <- r + (1 - 2*b1)*eps*gr(theta)
    theta <- theta + (0.5 - a1)*eps*r
    r <- r + b1*eps*gr(theta)
    theta <- theta + a1*eps*r
    
    ## verify valid trajectory
    H <- .calculate.H(theta=theta, r=r, fn=fn)
    s <- H-log(u) + delta.max > 0
    if(is.na(s) | is.nan(s)) s <- 0
    n <- log(u) <= H
    ## ## Useful code for debugging. Returns entire path to global env.
    ## if(!exists('theta.trajectory'))
    ##     theta.trajectory <<- theta
    ## else
    ##     theta.trajectory <<- rbind(theta.trajectory, theta)
    temp <- .calculate.H(theta=theta, r=r, fn=fn)-
      .calculate.H(theta=theta0, r=r0, fn=fn)
    alpha <- min(exp(temp),1)
    info$n.calls <- info$n.calls + 5
    return(list(theta.minus=theta, theta.plus=theta, theta.prime=theta, r.minus=r,
                r.plus=r, s=s, n=n, alpha=alpha, nalpha=1))
  } else {
    ## recursion - build left and right subtrees
    xx <- .buildtree.ths(theta=theta, r=r, u=u, v=v, j=j-1, eps=eps,
                     theta0=theta0, r0=r0, fn=fn, gr=gr, info=info)
    theta.minus <- xx$theta.minus
    theta.plus <- xx$theta.plus
    theta.prime <- xx$theta.prime
    r.minus <- xx$r.minus
    r.plus <- xx$r.plus
    alpha <- xx$alpha
    nalpha <- xx$nalpha
    s <- xx$s
    if(is.na(s) | is.nan(s)) s <- 0
    nprime <- xx$n
    ## If it didn't fail, update the above quantities
    if(s==1){
      if(v== -1){
        yy <- .buildtree.ths(theta=theta.minus, r=r.minus, u=u, v=v,
                         j=j-1, eps=eps, theta0=theta0, r0=r0,
                         fn=fn, gr=gr, info=info)
        theta.minus <- yy$theta.minus
        r.minus <- yy$r.minus
      } else {
        yy <- .buildtree.ths(theta=theta.plus, r=r.plus, u=u, v=v,
                         j=j-1, eps=eps, theta0=theta0, r0=r0,
                         fn=fn, gr=gr, info=info)
        theta.plus <- yy$theta.plus
        r.plus <- yy$r.plus
      }
      ## This isn't in the paper but if both slice variables failed,
      ## then you get 0/0. So I skip this test. Likewise if model
      ## throwing errors, don't keep that theta.
      nprime <- yy$n+ xx$n
      if(!is.finite(nprime)) nprime <- 0
      if(nprime!=0){
        ## choose whether to keep this theta
        if(runif(n=1, min=0, max=1) <= yy$n/nprime)
          theta.prime <- yy$theta.prime
      }
      ## check for valid proposal
      test <- .test.nuts(theta.plus=theta.plus,
                         theta.minus=theta.minus, r.plus=r.plus,
                         r.minus=r.minus)
      ## if(!test) warning(paste("U turn at j=", j))
      ## check if any of the stopping conditions were met
      s <- xx$s*yy$s*test
    }
    return(list(theta.minus=theta.minus, theta.plus=theta.plus,
                theta.prime=theta.prime,
                r.minus=r.minus, r.plus=r.plus, s=s, n=nprime,
                alpha=alpha, nalpha=1))
  }
}

#Finding epsilon
.find.epsilon.ths <- function(theta,  fn, gr, eps=1, verbose=TRUE){
  r <- rnorm(n=length(theta), mean=0, sd=1)
  ## Do one THREE STAGE step
  a1 <- 0.11888
  b1 <- 0.29620
  theta.new <- theta + a1*eps*r
  r.new <- r + b1*eps*gr(theta.new)
  theta.new <- theta.new + (0.5 - a1)*eps*r.new
  r.new <- r.new + (1 - 2*b1)*eps*gr(theta.new)
  theta.new <- theta.new + (0.5 - a1)*eps*r.new
  r.new <- r.new + b1*eps*gr(theta.new)
  theta.new <- theta.new + a1*eps*r.new
  
  H1 <- .calculate.H(theta=theta, r=r, fn=fn)
  H2 <- .calculate.H(theta=theta.new, r=r.new, fn=fn)
  a <- 2*(exp(H2)/exp(H1)>.5)-1
  ## If jumped into bad region, a can be NaN so setup algorithm to keep
  ## halving eps instead of throwing error
  if(!is.finite(a)) a <- -1
  k <- 1
  ## Similarly, keep going if there are infinite values
  while (!is.finite(H1) | !is.finite(H2) | a*H2-a*H1 > -a*log(2)) {
    eps <- (2^a)*eps
    ## Do one leapfrog step
    a1 <- 0.11888
    b1 <- 0.29620
    theta.new <- theta + a1*eps*r
    r.new <- r + b1*eps*gr(theta.new)
    theta.new <- theta.new + (0.5 - a1)*eps*r.new
    r.new <- r.new + (1 - 2*b1)*eps*gr(theta.new)
    theta.new <- theta.new + (0.5 - a1)*eps*r.new
    r.new <- r.new + b1*eps*gr(theta.new)
    theta.new <- theta.new + a1*eps*r.new
    H2 <- .calculate.H(theta=theta.new, r=r.new, fn=fn)
    k <- k+1
    if(k>50) {
      stop("More than 50 iterations to find reasonable eps. Model is likely misspecified or some other issue.")
    }
  }
  if(verbose) message(paste("Reasonable epsilon=", eps, "found after", k, "steps"))
  return(invisible(eps))
}

