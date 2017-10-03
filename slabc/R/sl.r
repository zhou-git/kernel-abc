## Synthetic likelihood approach to fitting non-linear dynamic models

sl.acf <- function(x,max.lag,coef) {
## `x' is a matrix containing replicate simulations in its columns.
## sl.acf turns these into acf's

  NAcode <- -1e70
  x[is.na(x)] <- NAcode
  
  acf <- matrix(0,max.lag+1,ncol(x))
  oo <- .C("slacf",acf=as.double(acf),x=as.double(x),as.integer(nrow(x)),as.integer(ncol(x)),
            as.integer(max.lag),as.double(NAcode),correlation=as.integer(coef),PACKAGE="slabc")

  acf <- matrix(oo$acf,max.lag+1,ncol(x))
  acf[acf == NAcode] <- NA
  acf

} ## end of sl.acf

nlar <- function(x,lag,power) {
## relatively efficient polynomial autoregression for multiple reps.
## each column of `x' is a replicate. 
## `lag[i]' is the lag for term i on rhs of autoregression
## `power[i]' is the power for term i on rhs of autoregression 
  beta <- matrix(0,length(lag),ncol(x))
  
  NAcode <- -1e70
  x[is.na(x)] <- NAcode  

  oo <- .C("slnlar",beta = as.double(beta), x = as.double(x),
            n=as.integer(nrow(x)),n.reps=as.integer(ncol(x)),n.terms=as.integer(length(lag)),
            as.integer(lag),as.integer(power),as.double(NAcode),PACKAGE="slabc")
  
   beta <- matrix(oo$beta,length(lag),ncol(x))

   beta
} ## end of nlar


order.dist <- function(x,z,np=3,diff=1) {
## Routine to obtain coefficients summarizing distribution of (differenced) columns
## of x, by regression of sorted differenced columns of x on sorted differenced z's. 
## regression is with order `np' polynomial (no intercept as all centred). `diff'
## is order of differencing to apply.
  
  beta <- matrix(0,np,ncol(x))
  oo <- .C("order_reg",beta=as.double(beta), as.double(x),as.double(z),n=as.integer(nrow(x)),
            as.integer(ncol(x)),as.integer(np),as.integer(diff),PACKAGE="slabc")
  
  beta <- matrix(oo$beta,np,ncol(x))

} ## end of order.dist


#### Model code 

ricker <- function(theta,e,u,burn.in) {
## e is noise matrix. ncol(e) is number of reps, nrow(e) is length of sim 
## including burn.in period, theta is parameter vector
## u is a sequence of uniform deviates, used for generating observed
## from expected while avoiding Monte Carlo error.

  n.t <- nrow(e)-burn.in
  n.reps <- ncol(e)
  n <- matrix(0,n.t,n.reps)
  oo <- .C("ricker",n=as.double(n),as.double(theta),as.double(e),as.integer(burn.in),as.integer(n.t),
                    as.integer(n.reps),as.double(1),PACKAGE="slabc")
  N <- exp(oo$n)           ## N|e
#  seed <- get(".Random.seed", envir = .GlobalEnv)
#  set.seed(10)
#  y <- rpois(N,theta[3]*N) ## data y where E(y|e) = N|e
#  assign(".Random.seed", seed, envir = .GlobalEnv)
  y <- qpois(u,theta[3]*N) ## data y where E(y|e) = N|e
  matrix(y,n.t,n.reps)
}

ricker_abc <- function(phi, log_r, sig_e, e, u, burn.in)
{
  n.t<-nrow(e) - burn.in
  n.reps <- ncol(e)
  n <- matrix(0, n.t, n.reps)
  oo<- .C("ricker_abc", n=as.double(n), as.double(log_r), as.double(sig_e), as.double(e), as.integer(burn.in),as.integer(n.t),
          as.integer(n.reps), as.double(1),  PACKAGE="slabc")
  N <- exp(oo$n)
  #   N[N==0] <- NA
  #   temp <- matrix(N,nrow=50)
  #   Nm<-temp[,colSums(!is.na(temp))>5]
  #   Nm[is.na(Nm)] <- 0
  #   n.col <- ncol(Nm)
  #   Nt <- c(Nm)
  #  phim <- rep(phi,times=1,each=50)
  #   n.eff <- length(Nt)
  #   phit <- phim[1:n.eff]
  #   ut <- u[1:n.eff]
  y <- qpois(u,theta[3]*N)## data y where E(y|e) = N|e
  matrix(y, n.t, n.reps)
}

ng.bf <- function(theta,lu,lu1,burn.in) {
## lu,lu1 are noise matrices of the same dimension.
## they contain logit(uniform) random deviates,
## fixed from rep to rep to avoid monte carlo error. 
## ncol(lu) is number of reps, nrow(lu) is length of sim 
## including burn.in period, theta is parameter vector

## first transform the logit(uniform) deviates in lu and lu1 to
## gamma deviates with mean 1 and given variance, and do it quickly!
  var <- theta[4];var1 <- theta[6]

  luk <- seq(min(lu),max(lu),length=200) ## even spacing on logit(prob) scale
  uk <- exp(luk);uk <- uk/(1+uk)          ## on prob scale
  gk <- qgamma(uk,shape=1/var,scale=var) ## reference quantiles
  gk[!is.finite(gk)] <- 0                ## can underflow badly at lower end
  e <- approx(x=luk,y=gk,xout=lu)$y      ## fast approximate version of qgamma(u,...) 

  luk <- seq(min(lu1),max(lu1),length=200)
  uk <- exp(luk);uk <- uk/(1+uk)
  gk <- qgamma(uk,shape=1/var1,scale=var1)
  gk[!is.finite(gk)] <- 0 
  e1 <- approx(x=luk,y=gk,xout=lu1)$y 

  n.t <- nrow(lu)-burn.in
  n.reps <- ncol(lu)
  n <- matrix(0,n.t,n.reps)
  oo <- .C("ng_bf",n=as.double(n),as.double(theta),as.double(e),as.double(e1),as.integer(burn.in),as.integer(n.t),
                    as.integer(n.reps),PACKAGE="slabc")
  matrix(oo$n,n.t,n.reps)
}

des.bf <- function(theta,burn.in=100,n.t=100,n.rep=500) {
## Blowfly model with demographic stochasticity + environmental stochasticity
## as required by a referee. No real reason to use this rather than 
## ng.bf: demographic stochasticity is trivial here. Each column is a replicate.
## theta <- c(.16,6.5,400,14) 
## pname <- names(theta) <- c("delta","P","N0","tau")
## N <- ds.bf(theta) 

  delta <- theta[1]
  P <- theta[2]
  N0 <- theta[3]
  tau <- theta[4]
  var.P <- theta[5];
  var.d <- theta[6]
  ## gamma deviates for env stoch.... 
  gp <- matrix(rgamma((n.t+burn.in)*n.rep,shape=1/var.P,scale=var.P),burn.in+n.t,n.rep)
  gd <- matrix(rgamma((n.t+burn.in)*n.rep,shape=1/var.d,scale=var.d),burn.in+n.t,n.rep)
  N <- matrix(1:n.rep*2+500,burn.in+n.t,n.rep)
  for (i in 1:tau) { 
    B <- P * exp(-N[1,]/N0)*gp[i,]; ## birth rate
    D <- exp(-delta*gd[i,]) ## survival rate
    N[i+1,] <- rpois(N[1,],N[1,]*B) + rbinom(N[i,],N[i,],D) ## Birth + Survival    
  }
  for (i in (tau+1):(burn.in-1)) {
    B <- P * exp(-N[i-tau,]/N0)*gp[i,]; ## birth rate
    D <- exp(-delta*gd[i,]) ## survival rate
    N[i+1,] <- rpois(N[1,],N[i-tau,]*B) + rbinom(N[i,],N[i,],D) ## Birth + Survival    
  }
  for (i in burn.in:(burn.in+n.t-1)) {
    B <- P * exp(-N[i-tau,]/N0)*gp[i,]; ## birth rate
    D <- exp(-delta*gd[i,]) ## survival rate
    N[i+1,] <- rpois(N[1,],N[i-tau,]*B) + rbinom(N[i,],N[i,],D) ## Birth + Survival    
  }
  N[(burn.in+1):(n.t+burn.in),]
}


ds.bf <- function(theta,burn.in=100,n.t=100,n.rep=500) {
## Blowfly model with pure demographic stochasticity. Each column is a replicate.
## theta <- c(.16,6.5,400,14) 
## pname <- names(theta) <- c("delta","P","N0","tau")
## N <- ds.bf(theta) 

  delta <- theta[1]
  P <- theta[2]
  N0 <- theta[3]
  tau <- theta[4]
  N <- matrix(1:n.rep*2+500,burn.in+n.t,n.rep)
  D <- exp(-delta)
  for (i in 1:tau) { 
    B <- P * exp(-N[1,]/N0); ## birth rate
    N[i+1,] <- rpois(N[1,],N[1,]*B) + rbinom(N[i,],N[i,],D) ## Birth + Survival    
  }
  for (i in (tau+1):(burn.in-1)) {
    B <- P * exp(-N[i-tau,]/N0); ## birth rate
    N[i+1,] <- rpois(N[1,],N[i-tau,]*B) + rbinom(N[i,],N[i,],D) ## Birth + Survival    
  }
  for (i in burn.in:(burn.in+n.t-1)) {
    B <- P * exp(-N[i-tau,]/N0); ## birth rate
    N[i+1,] <- rpois(N[1,],N[i-tau,]*B) + rbinom(N[i,],N[i,],D) ## Birth + Survival    
  }
  N[(burn.in+1):(n.t+burn.in),]
}


bup.para <- function(theta,e,u=0,burn.in) {
## e is noise matrix. ncol(e) is number of reps, nrow(e) is length of sim 
## including burn.in period, theta is parameter vector
## u is a sequence of uniform deviates, used for generating observed
## from expected while avoiding Monte Carlo error.

  n.t <- nrow(e)-burn.in
  n.reps <- ncol(e)
  n <- matrix(0,n.t,n.reps)
  oo <- .C("bup_par",n=as.double(n),as.double(theta),as.double(e),as.integer(burn.in),as.integer(n.t),
                    as.integer(n.reps),PACKAGE="slabc")
  
  matrix(oo$n^.25+qnorm(u,0)*theta[6],n.t,n.reps) ## transformed response with obs error
}



## preliminary code

print.sl.version <- function()
{ library(help=slabc)$info[[1]] -> version
  version <- version[pmatch("Version",version)]
  um <- strsplit(version," ")[[1]]
  version <- um[nchar(um)>0][2]
  cat(paste("This is sl ",version,"\n"),sep="")
}


.onAttach <- function(...) { 
  print.sl.version()
}

.onUnload <- function(libpath) library.dynam.unload("slabc", libpath)

.onLoad <- function(lib,pkg) {
   library.dynam("slabc", pkg, lib)
}

.First.lib <- function(lib, pkg) {
  library.dynam("slabc", pkg, lib)
  print.sl.version()
}
