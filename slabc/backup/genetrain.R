library(slabc)
library(parallel)
set.seed(3);
 cl <- makeForkCluster(16 )
snum<-1e+06
##Parameter setting 
n.total <- snum*2
burn.in <- 50
n.t     <- 50

train.log_r     <- runif(n.total, min = 0, max = 10)
train.log_sig_e <- runif(n.total, min = log(0.1), max = 0)
train.sig_e     <- exp(train.log_sig_e)
train.phi       <- runif(n.total, min = 0, max <- 30)

e <- matrix(rnorm((burn.in+n.t) * n.total), burn.in + n.t, n.total)
u <- runif(n.t * n.total) 
n.total <- snum
## simulate data ....
train.yt <- ricker_abc(train.phi, train.log_r, train.sig_e, e, u, burn.in)
train.yt[train.yt==0]<-NA

ttable <-colSums(!is.na(train.yt))>5
train.ytm <- train.yt[,ttable]
train.ytm[is.na(train.ytm)] <- 0

train.log_rf <- train.log_r[ttable]
train.log_sig_ef<- train.log_sig_e[ttable]
train.sig_ef<- train.sig_e[ttable]
train.phif<- train.phi[ttable]

train.y<-train.ytm[,1:n.total]
train.log_r <- train.log_rf[1:n.total]
train.log_sig_e<- train.log_sig_ef[1:n.total]
train.sig_e<- train.sig_ef[1:n.total]
train.phi<- train.phif[1:n.total]

## Now assemble the relevant statistics
acf.train.y  <- sl.acf(train.y,max.lag=5,0)
accf.train.y <- sl.acf(train.y,max.lag=5,1)
b0.train.y <- nlar(train.y^.3,lag=c(1,1),power=c(1,2))
#b1.train.y <- order.dist(train.y,train.y,np=3,diff=1)  
colVars <- function(x, na.rm=FALSE, dims=1, unbiased=TRUE, SumSquares=FALSE,
                    twopass=FALSE) {
  if (SumSquares) return(colSums(x^2, na.rm, dims))
  N <- colSums(!is.na(x), FALSE, dims)
  Nm1 <- if (unbiased) N-1 else N
  if (twopass) {x <- if (dims==length(dim(x))) x - mean(x, na.rm=na.rm) else
    sweep(x, (dims+1):length(dim(x)), colMeans(x,na.rm,dims))}
  (colSums(x^2, na.rm, dims) - colSums(x, na.rm, dims)^2/N) / Nm1
}


train.yo <- apply(train.y, 2, sort)
train.yad <- train.y[2:n.t, ] - train.y[1:n.t-1, ]  #time difference
train.yado <- apply(train.yad, 2, sort)

## combine the statistics...
train.sy.E0 <- rbind(acf.train.y,
               b0.train.y, 
               colMeans(train.y),
               colSums(train.y==0)
)

train.sy.E1 <- rbind(acf.train.y, b0.train.y, 
               colMeans(train.y), colSums(train.y==0),
               colSums(train.y==1), colSums(train.y==2),
               colSums(train.y==3), colSums(train.y==4),
               log(colMeans(train.y)), colVars(train.y),
               log(colSums(train.y^2)),
               log(colSums(train.y^3)),
               log(colSums(train.y^4)),
               log(colSums(train.y^5)),
               log(colSums(train.y^6)),
               accf.train.y[2:6,]
)
train.sy.E2 <- rbind(acf.train.y, b0.train.y, 
               colMeans(train.y), colSums(train.y==0),
               colSums(train.y==1), colSums(train.y==2),
               colSums(train.y==3), colSums(train.y==4),
               log(colMeans(train.y)), colVars(train.y),
               log(colSums(train.y^2)),
               log(colSums(train.y^3)),
               log(colSums(train.y^4)),
               log(colSums(train.y^5)),
               log(colSums(train.y^6)),
               accf.train.y[2:6,], train.y, train.yo, train.y^2, train.yo^2,
               log(1 + train.y), log(1 + train.yo),
               train.yad^2, train.yado^2
)

save(train.log_r, train.log_sig_e, train.phi, train.y, train.sy.E0, train.sy.E1, train.sy.E2, file = "ricketrain10e06.Rdata")
