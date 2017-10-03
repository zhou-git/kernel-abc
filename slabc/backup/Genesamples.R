library(slabc)
# library(parallel)
set.seed(3);
# cl <- makeForkCluster(16 )
snum<-1e+06
##Parameter setting 
n.total <- snum*2

burn.in <- 50
n.t     <- 50

log_r     <- runif(n.total, min = 0, max = 10)
log_sig_e <- runif(n.total, min = log(0.1), max = 0)
sig_e     <- exp(log_sig_e)
phi       <- runif(n.total, min = 0, max <- 30)

e <- matrix(rnorm((burn.in+n.t) * n.total), burn.in + n.t, n.total)
u <- runif(n.t * n.total) 
n.total <- snum
## simulate data ....
Yt <- ricker_abc(phi, log_r, sig_e, e, u, burn.in)
Yt[Yt==0]<-NA

ttable <-colSums(!is.na(Yt))>5
Ytm <- Yt[,ttable]
Ytm[is.na(Ytm)] <- 0

log_rf <- log_r[ttable]
log_sig_ef<- log_sig_e[ttable]
sig_ef<- sig_e[ttable]
phif<- phi[ttable]

Y<-Ytm[,1:n.total]
log_r <- log_rf[1:n.total]
log_sig_e<- log_sig_ef[1:n.total]
sig_e<- sig_ef[1:n.total]
phi<- phif[1:n.total]
## Now assemble the relevant statistics
acf.Y  <- sl.acf(Y,max.lag=5,0)
accf.Y <- sl.acf(Y,max.lag=5,1)
b0.Y <- nlar(Y^.3,lag=c(1,1),power=c(1,2))
#b1.Y <- order.dist(Y,y,np=3,diff=1)  
colVars <- function(x, na.rm=FALSE, dims=1, unbiased=TRUE, SumSquares=FALSE,
                    twopass=FALSE) {
  if (SumSquares) return(colSums(x^2, na.rm, dims))
  N <- colSums(!is.na(x), FALSE, dims)
  Nm1 <- if (unbiased) N-1 else N
  if (twopass) {x <- if (dims==length(dim(x))) x - mean(x, na.rm=na.rm) else
    sweep(x, (dims+1):length(dim(x)), colMeans(x,na.rm,dims))}
  (colSums(x^2, na.rm, dims) - colSums(x, na.rm, dims)^2/N) / Nm1
}


Yo <- apply(Y, 2, sort)
Yad <- Y[2:n.t, ] - Y[1:n.t-1, ]  #time difference
Yado <- apply(Yad, 2, sort)

## combine the statistics...
sY.E0 <- rbind(acf.Y,
               b0.Y, 
               colMeans(Y),
               colSums(Y==0)
)

sY.E1 <- rbind(acf.Y, b0.Y, 
               colMeans(Y), colSums(Y==0),
               colSums(Y==1), colSums(Y==2),
               colSums(Y==3), colSums(Y==4),
               log(colMeans(Y)), colVars(Y),
               log(colSums(Y^2)),
               log(colSums(Y^3)),
               log(colSums(Y^4)),
               log(colSums(Y^5)),
               log(colSums(Y^6)),
               accf.Y[2:6,]
)
sY.E2 <- rbind(acf.Y, b0.Y, 
               colMeans(Y), colSums(Y==0),
               colSums(Y==1), colSums(Y==2),
               colSums(Y==3), colSums(Y==4),
               log(colMeans(Y)), colVars(Y),
               log(colSums(Y^2)),
               log(colSums(Y^3)),
               log(colSums(Y^4)),
               log(colSums(Y^5)),
               log(colSums(Y^6)),
               accf.Y[2:6,], Y, Yo, Y^2, Yo^2,
               log(1 + Y), log(1 + Yo),
               Yad^2, Yado^2
)

save(log_r, log_sig_e, phi, Y, sY.E0, sY.E1, sY.E2, file = "rickesample.Rdata")

