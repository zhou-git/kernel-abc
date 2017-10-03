## routine's for transformation of statistics to better meet normality assumptions,
## and for checking the MVN approximation


get.trans <- function(S) {
## S is an ns by n.reps matrix of statistics. This routine works through
## its rows finding piecewise linear transformations to normality, by 
## interactive use of `locator'.
  op <- par(mfrow=c(1,1))
  if (!is.matrix(S)) S <- matrix(S,1,length(S))
  ns <- nrow(S)    ## the number of statistics
  n.rep <- ncol(S) ## the number of replicates
  z <- qnorm((1:n.rep-.5)/n.rep) ## standard normal quantiles
  trans <- list()
  for (i in 1:ns) { ## loop through stats...
    plot(sort(S[i,]),z)
    tr <- locator(,type="l",col=2)
    if (length(tr$x)<2) { 
      warning("no transform");
      trans[[i]] <- NULL  
    } else
    ## now extend end segments
    if (length(tr$x)==2) { ## single segment --- extend both ends
      xr <- tr$x[2]-tr$x[1]
      slope <- (tr$y[2]-tr$y[1])/xr
      tr$x[1] <- tr$x[1] - 1000*xr 
      tr$y[1] <- tr$y[1] - slope*1000*xr
      tr$x[2] <- tr$x[2] + 1000*xr
      tr$y[2] <- tr$y[2] + slope*1000*xr
      trans[[i]] <- tr      
    } else { ## extend end segments
      xr <- max(tr$x) - min(tr$x)
      slope <- (tr$y[2]-tr$y[1])/(tr$x[2]-tr$x[1])
      tr$x[1] <- tr$x[1] - 1000*xr 
      tr$y[1] <- tr$y[1] - slope*1000*xr
      nt <- length(tr$x)
      slope <- (tr$y[nt]-tr$y[nt-1])/(tr$x[nt]-tr$x[nt-1])
      tr$x[nt] <- tr$x[nt] + 1000*xr
      tr$y[nt] <- tr$y[nt] + slope*1000*xr
      trans[[i]] <- tr     
    }
  }
  trans[[ns+1]] <- NA
  par(op)
  trans ## the transformation object
}

trans.stat <- function(S,trans) {
## apply a piecewise linear `trans' object to the rows 
## of a statistics matrix, S
  if (!is.matrix(S)) S <- matrix(S,length(S),1)
  for (i in 1:nrow(S)) {
    if (!is.null(trans[[i]]))
    S[i,] <- approx(trans[[i]]$x,trans[[i]]$y,S[i,],rule=2)$y
  }
  if (ncol(S)==1) S <- as.numeric(S)
  S
}


not.sq <- function(x,alpha= .5,d0 = 5) {
## not.sq(x) = x^2           x <= d0
##             k*x^alpha + c otherwise
## k and c are chosen to give continuity to first derivative
  ind <- x < d0
  x[ind] <- x[ind]^2
  k <- 2*d0^(2-alpha)/alpha
  cc <- d0^2 - k * d0^alpha
  x[!ind] <- k * x[!ind]^ alpha + cc  
  x
}


MVN.check <- function(S,s=NULL,cex.axis=1,cex.lab=1) {
## graphical check of multivariate normality, from
## Krzanowski (1988) 7.5
  p <- nrow(S)
  n <- ncol(S)
  if (n<10*p) warning("You don't really have enough reps for this approach")
  
  ps <- s
  for (i in 1:nrow(S)) ps[i] <- sum(S[i,]<s[i])/ncol(S)

  um <- robust.vcov(S)
  ms <- as.numeric(um$mY)
  S <- S-ms
  ## Malahanobis for each column of S
  z <- colSums(S*(t(um$E)%*%um$E%*%S))
  
  q <- log(qchisq((1:n-.5)/n,df=p))
  z <- log(sort(z))
  
  plot(q,z,type="l",col="grey",
       xlab="log theoretical quantiles",ylab="log observed quantiles",
       cex.axis=cex.axis,cex.lab=cex.lab)
  points(q,z,pch=".")
  abline(0,1,col=2)
  cat("\nproportion |log(z)-log(q)|>.25 = ",sum(abs(z-q)>.25)/n,"\n")
  
  if (!is.null(s)) { ## QQ plot for observed stats
    z <- um$E%*%(s-ms)
    q <- sum(z^2)
    abline(h=log(q),lty=2)
  }


  ## Marginal stats

  for (i in 1:nrow(S)) S[i,] <- S[i,]/um$sd[i]
  n <- ncol(S)
  z <- qnorm((1:n-.5)/n)
  rz <- range(z)
  rz <- c(rz[1]-2,rz[2]+2)
  
  plot(z,sort(S[1,]),type="l",col="grey",ylim=rz,
       xlab="N(0,1) quantiles",ylab="marginal quantiles",
       cex.axis=cex.axis,cex.lab=cex.lab)
  points(z,sort(S[1,]),pch=".")
  for (i in 2:nrow(S)) { 
    lines(z,sort(S[i,]),col="grey")
    points(z,sort(S[i,]),pch=".")
  }
  abline(0,1,col=2)  

  if (!is.null(s)) { ## QQ plot for observed stats
    z <- um$E%*%(s-ms)
    qqnorm(z,cex.axis=cex.axis,cex.lab=cex.lab);qqline(z,col=2)
  }

  ps
 
}