colVars <- function(x, na.rm=FALSE, dims=1, unbiased=TRUE, SumSquares=FALSE,
                    twopass=FALSE) {
  if (SumSquares) return(colSums(x^2, na.rm, dims))
  N <- colSums(!is.na(x), FALSE, dims)
  Nm1 <- if (unbiased) N-1 else N
  if (twopass) {x <- if (dims==length(dim(x))) x - mean(x, na.rm=na.rm) else
    sweep(x, (dims+1):length(dim(x)), colMeans(x,na.rm,dims))}
  (colSums(x^2, na.rm, dims) - colSums(x, na.rm, dims)^2/N) / Nm1
}

ricker_rej <- function(i, n.t, n.acc, n.total, n.obs, c.E, trans=NULL, stats=FALSE) {
  ##rejection method for ricker abc.
  
  ## simulate from model with current parameter values
  #load("rickerObs100.Rdata")
  #load("rickerSample10000.Rdata")
  
#   if (!is.matrix(obs.y)) {
#     obs.y <- matrix(obs.y, n.t, n.obs)
#   }
#   y <- matrix(obs.y[, i], length(obs.y[, 1]), 1)
  y<-obs.y
  obs.log_r.i <-obs.log_r[i]
  obs.log_sig_e.i <- obs.log_sig_e[i]
  obs.phi.i <- obs.phi[i]
  theta.obs <- c(obs.log_r.i, obs.log_sig_e.i, obs.phi.i)
  
  ## Now assemble the relevant statistics
  #acf.Y  <- sl.acf(Y,max.lag=5,0)
  acf.y  <- sl.acf(y,max.lag=5,0)
#  accf.Y <- sl.acf(Y,max.lag=5,1)
  accf.y <- sl.acf(y,max.lag=5,1)
  
 # b0.Y <- nlar(Y^.3,lag=c(1,1),power=c(1,2))
  b0.y <- nlar(y^.3,lag=c(1,1),power=c(1,2))
  
#   b1.Y <- order.dist(Y,y,np=3,diff=1)
  b1.y <- order.dist(y,y,np=3,diff=1)   
  
#  Yo <- apply(Y, 2, sort)
  yo <- apply(y, 2, sort)
 # Yad <- Y[2:n.t, ] - Y[1:n.t-1, ]  #time difference
  yad <- y[2:n.t, ] - y[1:n.t-1, ]
 # Yado <- apply(Yad, 2, sort)
  yado <- apply(yado, 2, sort)
  
  ## combine the statistics...
#   sY.E0 <- rbind(sY.E0,b1.Y)
#   sY.E1 <- rbind(sY.E1,b1.Y)
#   sY.E2 <- rbind(sY.E2,b1.Y)

#   sy.E0 <- c(as.numeric(acf.y),
#              as.numeric(b0.y),
#              mean(y),sum(y==0),
#              as.numeric(b1.y)
#   )
  sy.E0 <- rbind(acf.y,
              b0.y,
              colMeans(y),
              colSums(y==0), 
              b1.y
  )
#   sy.E1 <- c(as.numeric(acf.y),
#              as.numeric(b0.y),
#              mean(y), sum(y==0),
#              sum(y==1), sum(y==2),
#              sum(y==3), sum(y==4),
#              as.numeric(log(mean(y))),
#              as.numeric(var(y)),
#              as.numeric(log(sum(y^2))),
#              as.numeric(log(sum(y^3))),
#              as.numeric(log(sum(y^4))),
#              as.numeric(log(sum(y^5))),
#              as.numeric(log(sum(y^6))),
#              as.numeric(accf.y),
#              as.numeric(b1.y)
#     )
  sy.E1 <- rbind(acf.y, b0.y, 
             colMeans(y), colSums(y==0),
             colSums(y==1), colSums(y==2),
             colSums(y==3), colSums(y==4),
             log(colMeans(y)), colVars(y),
             log(colSums(y^2)),
             log(colSums(y^3)),
             log(colSums(y^4)),
             log(colSums(y^5)),
             log(colSums(y^6)),
             accf.y, b1.y
    )
#   sy.E2 <- c(as.numeric(acf.y),
#              as.numeric(b0.y),
#              mean(y), sum(y==0),
#              sum(y==1), sum(y==2),
#              sum(y==3), sum(y==4),
#              as.numeric(log(mean(y))),
#              as.numeric(var(y)),
#              as.numeric(log(sum(y^2))),
#              as.numeric(log(sum(y^3))),
#              as.numeric(log(sum(y^4))),
#              as.numeric(log(sum(y^5))),
#              as.numeric(log(sum(y^6))),
#              as.numeric(accf.y),
#              as.numeric(y),
#              as.numeric(yo),
#              as.numeric(y^2),
#              as.numeric(yo^2),
#              as.numeric(log(1+y)),
#              as.numeric(log(1+yo)),
#              as.numeric(yad^2),
#              as.numeric(yado^2),
#              as.numeric(b1.y)
#     )
  sy.E2 <- rbind(acf.y, b0.y, 
                 colMeans(y), colSums(y==0),
                 colSums(y==1), colSums(y==2),
                 colSums(y==3), colSums(y==4),
                 log(colMeans(y)), colVars(y),
                 log(colSums(y^2)),
                 log(colSums(y^3)),
                 log(colSums(y^4)),
                 log(colSums(y^5)),
                 log(colSums(y^6)),
                 accf.y, y, yo, y^2, yo^2,
                 log(1 + y), log(1 + yo),
                 yad^2, yado^2, b1.y
    )
  
  if (!is.null(trans)) {
    sy.E0 <- trans.stat(sy.E0,trans)
    sY.E0 <- trans.stat(sY.E0,trans)
  }
  
  ## get the log synthetic likelihood
  if(c.E == 1) {
    sY.EC <- sY.E0
    sy.EC <- sy.E0
  }else if(c.E == 2) {
    sY.EC <- sY.E1
    sy.EC <- sy.E1
  }else {
    sY.EC <- sY.E2
    sy.EC <- sy.E2
  }
  
  sY.EC <- sY.EC[,is.finite(colSums(sY.EC))]
  if (stats) {  ## return statistics
    attr(sY.EC,"observed") <- sy.EC
    return(sY.EC) 
  }
  
  #er <- robust.vcov(sY.EC)
  Ymad <- apply(sY.EC, 1, mad)
  Ymad[which(Ymad == 0)] <- 1
  minusSS <- apply(sY.EC, 2, '-', sy.EC)/Ymad
  rss <- sqrt(colSums((minusSS^2)/length(sy.EC)))
  index.total <- sort(rss, index.return = TRUE)
  index.acc <- index.total$ix[1:100]
  
  #output
  log_r.mean <- mean(log_r[index.acc])
  log_e.mean <- mean(log_sig_e[index.acc])
  phi.mean <- mean(phi[index.acc])
  theta.mean <- c(log_r.mean, log_e.mean, phi.mean)
  rsse <- dist(rbind(theta.obs, theta.mean))
  
  #rss <- sum((er$E%*%(sy.EC-er$mY))^2)
  #ll <- -rss/2 - er$half.ldet.V
  
}