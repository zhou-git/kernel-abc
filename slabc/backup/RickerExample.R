library(sl)

set.seed(3)
n.t <- 50 ## length of data run
burn.in <- 50 ## burnin length for model
n.reps <- 500 ## slow mixing if this lowered 

## simulate random deviates to re-use for fitting...
e <- matrix(rnorm((burn.in+n.t)*n.reps),burn.in+n.t,n.reps)
u <- runif(n.t*n.reps) 

## simulate data ....
theta <- c(3.8,0.3,10)
names(theta) <- c("r","sig.e","phi")
Y <- ricker(theta,e,u,burn.in)
y <- Y[,1]

## plot example replicates from model... 
par(mfrow=c(2,2))
plot(1:n.t,Y[,1],xlab="t",ylab="y",type="l")
for (i in 2:5) lines(1:n.t,Y[,i],col=i)

## Check normality of statistics...
sY <- ricker.ll(theta,y,e,u,burn.in,stats=TRUE)
sy <- attr(sY,"observed")
par(mfrow=c(2,2))
MVN.check(sY,sy)


## Transects of the log synthetic likelihood...

par(mfrow=c(2,2))
lo <- c(3,.1,4); hi <- c(5,.8,20)
for (k in 1:3) {
  ll <- x <- seq(lo[k],hi[k],length=100)
  thetap <- theta
  for (i in 1:100) {
    thetap[k] <- x[i]
    ll[i] <- ricker.ll(thetap,y,e,u,burn.in)
  }
  plot(x,ll,type="l")
}

########################################################
## The MCMC chain exploring the synthetic log likelihood
########################################################

rktrans <- NULL ## default value (stats used raw)

## load transformations of statistics obtained after a pilot run with 
## untransformed statistics. Comment out this line to use 
## statistics untransformed....

data(rktrans) 

set.seed(3) ## comment out for replicates!!

n.mc <- 600 ## value for package checking
## Not run: 
n.mc <- 50000 ## length of MCMC run

## End(Not run)
theta <- c(4,0.1,6) ## gets dynamics somewhere in right ball park
n.para <- length(theta)
th <- matrix(0,n.para,n.mc) ## storage for chain results
llr <- rep(NA,n.mc)         ## storage for log synthetic likelihood
th[,1] <- log(theta)        ## initial state
prop.sd <- c(.02,.1,.05)    ## proposal standard deviations


llr[1] <- ricker.ll(exp(th[,1]),y,e,u,burn.in=burn.in,trans=rktrans)
reject <- 0
uni <- runif(n.mc)
pname <- c("r","sig.e","phi")

for (i in 2:n.mc) { ## MCMC loop
  th[,i] <- th[,i-1] + rnorm(n.para)*prop.sd*2
  
  llr[i] <- ricker.ll(exp(th[,i]),y,e,u,burn.in,trans=rktrans)
  
  alpha <- min(1,exp(llr[i]-llr[i-1]))
  
  if (uni[i]>alpha) { ## reject
    th[,i] <- th[,i-1]
    llr[i] <- llr[i-1]
    reject <- reject + 1
  }
  
  if (i%%200==0) { ## plot the chains so far (and report acceptance rate)
    par(mfrow=c(2,2))
    start <- 1
    for (j in 1:n.para) { 
      xx <- exp(th[j,1:i])
      if (i < 5000) ylim=range(xx) else ylim=range(xx[(i/2):i])
      plot(1:i,xx,type="l",ylab=pname[j],ylim=ylim)
    }
    xx <- llr[1:i]
    if (i < 5000) ylim=range(xx) else ylim=range(xx[(i/2):i])
    plot(1:i,xx,type="l",ylab="log l_s",ylim=ylim,main=round(1-reject/i,digits=2))
  }
} ## end of the MCMC chain

1-reject/n.mc ## acceptance rate


start <- round(n.mc/2)
for (j in 1:n.para) hist(exp(th[j,start:n.mc]),xlab=pname[j]); 

hist(llr[start:n.mc],xlab="log l_s")

## Do some "pure" likelihood based inference...

th <- rbind(th,llr)
rownames(th) <- c("r","sig","phi","ll")

ml <- chain2ll(th,start=start)

## normality check...
sY <- ricker.ll(exp(ml$mle),y,e,u,burn.in,trans=rktrans,stats=TRUE)
sy <- attr(sY,"observed")
MVN.check(sY,sy)


## Now examine the joint density of the data and random effects 
## as a function of r and random effects...
n.t <- 50 
burn.in <- 50
n.reps <- 1  
e <- matrix(rnorm((burn.in+n.t)*n.reps),burn.in+n.t,n.reps)
u <- runif(n.t*n.reps)
theta <- c(3.8,0.3,10)
names(theta) <- c("r","sig.e","phi")
y <- ricker(theta,e,u,burn.in)
ricker.fey(theta,y,e)
n <- 1000
l.fey <- r <- seq(2.5,4.5,length=n)
th <- theta
for (i in 1:n) { 
  th[1] <- r[i]
  l.fey[i] <- ricker.fey(th,y,e)
}
plot(r,l.fey,type="l")

fey <- e1 <- seq(-2,2,length=n)
th <- theta;ep <- e
for (i in 1:n) { 
  ep[1] <- e1[i]
  fey[i] <- ricker.fey(th,y,ep)
}
plot(e1,fey,type="l")
