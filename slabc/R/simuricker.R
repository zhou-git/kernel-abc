
simuRicker <- function (param) {
  if(length(param)>3)
  set.seed(param[1])##???
  ##Parameter setting
  log_r=param[2]
  log_sig_e=param[3]
  phi=param[4]
  ss_set=1
  snum=1
  n.total <- snum
  
  burn.in <- 50
  n.t     <- 50

  sig_e     <- exp(log_sig_e)
  
  e <- matrix(rnorm((burn.in+n.t) * n.total), burn.in + n.t, n.total)
  u <- runif(n.t * n.total) 
  n.total <- snum
  theta=c(log_r,sig_e,phi)
  ## simulate data ....
  Y <- ricker(theta, e, u, burn.in)


  # Yt[Yt==0]<-NA
  # 
  # ttable <-colSums(!is.na(Yt))>5
  # Ytm <- Yt[,ttable]
  # Ytm[is.na(Ytm)] <- 0
  # 
  # log_rf <- log_r[ttable]
  # log_sig_ef<- log_sig_e[ttable]
  # sig_ef<- sig_e[ttable]
  # phif<- phi[ttable]
  # 
  # Y<-Ytm[,1:n.total]
  # log_r <- log_rf[1:n.total]
  # log_sig_e<- log_sig_ef[1:n.total]
  # sig_e<- sig_ef[1:n.total]
  # phi<- phif[1:n.total]
  ## Now assemble the relevant statistics
  acf.Y  <- sl.acf(Y,max.lag=5,0)#dim:6
  accf.Y <- sl.acf(Y,max.lag=5,1)#6
  Yad = Y[2:n.t, ] - Y[1:n.t-1, ]
  b0.Y <- nlar(Y^.3,lag=c(1,1),power=c(1,2))#2
  b1.Y <- order.dist(Y,Yad,np=3,diff=1)  #3

  Yo <- apply(Y, 2, sort)
    #time difference
  Yad=matrix(Yad,nrow = 49, ncol = n.total)
  Yado <- apply(Yad, 2, sort)
  
  ## combine the statistics...
  if(ss_set==0){
    sY <- c(as.vector(acf.Y),
               as.vector(b0.Y), 
               as.vector(b1.Y),
                   colMeans(Y),
                   colSums(Y==0)
    )
  }
  if(ss_set==1){
    sY <- c(as.vector(acf.Y), as.vector(b0.Y), as.vector(b1.Y),
                   colMeans(Y), colSums(Y==0),
                   colSums(Y==1), colSums(Y==2),
                   colSums(Y==3), colSums(Y==4),
                   log(colMeans(Y)), 
                   log(colSums(Y^2)),
                   log(colSums(Y^3)),
                   log(colSums(Y^4)),
                   log(colSums(Y^5)),
                   log(colSums(Y^6)),
               as.vector(accf.Y[2:6,]))
  }
  if(ss_set==2){
    sY <- c(as.vector(acf.Y), as.vector(b0.Y),as.vector(b1.Y),
                   colMeans(Y), colSums(Y==0),
                   colSums(Y==1), colSums(Y==2),
                   colSums(Y==3), colSums(Y==4),
                   log(colMeans(Y)), 
                   log(colSums(Y^2)),
                   log(colSums(Y^3)),
                   log(colSums(Y^4)),
                   log(colSums(Y^5)),
                   log(colSums(Y^6)),
                   as.vector(accf.Y[2:6,]), as.vector(Y),as.vector(Yo), as.vector(Y^2), as.vector(Yo^2),
                   as.vector(log(1 + Y)), as.vector(log(1 + Yo)),
                   as.vector(Yad^2), as.vector(Yado^2)
    )
  }

  return(sY)
}

simRicker_gkdr <- function(snum,ss_set){
  
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
  train.yt <- ricker_abc(phi, log_r, sig_e, e, u, burn.in)
  train.yt[train.yt==0]<-NA
  
  ttable <-colSums(!is.na(train.yt))>5
  train.ytm <- train.yt[,ttable]
  train.ytm[is.na(train.ytm)] <- 0
  
  train.log_rf <- log_r[ttable]
  train.log_sig_ef<- log_sig_e[ttable]
  train.sig_ef<- sig_e[ttable]
  train.phif<- phi[ttable]
  
  train.y<-train.ytm[,1:n.total]
  train.log_r <- train.log_rf[1:n.total]
  train.log_sig_e<- train.log_sig_ef[1:n.total]
  train.sig_e<- train.sig_ef[1:n.total]
  train.phi<- train.phif[1:n.total]
  
  ## Now assemble the relevant statistics
  acf.train.y  <- sl.acf(train.y,max.lag=5,0)
  accf.train.y <- sl.acf(train.y,max.lag=5,1)
  b0.train.y <- nlar(train.y^.3,lag=c(1,1),power=c(1,2))
  b1.train.y <- order.dist(train.y,train.y,np=3,diff=1)  

  train.yo <- apply(train.y, 2, sort)
  train.yad <- train.y[2:n.t, ] - train.y[1:n.t-1, ]  #time difference
  train.yado <- apply(train.yad, 2, sort)
  
  ## combine the statistics...
  if(ss_set==0){
  E <- rbind(acf.train.y,
                       b0.train.y, b1.train.y,
                       colMeans(train.y),
                       colSums(train.y==0)
  )
  }
  if(ss_set==1){
  E <- rbind(acf.train.y, b0.train.y, b1.train.y,
                       colMeans(train.y), colSums(train.y==0),
                       colSums(train.y==1), colSums(train.y==2),
                       colSums(train.y==3), colSums(train.y==4),
                       log(colMeans(train.y)), 
                       log(colSums(train.y^2)),
                       log(colSums(train.y^3)),
                       log(colSums(train.y^4)),
                       log(colSums(train.y^5)),
                       log(colSums(train.y^6)),
                       accf.train.y[2:6,]
  )
  }
  if(ss_set==2){
  E <- rbind(acf.train.y, b0.train.y, b1.train.y,
                       colMeans(train.y), colSums(train.y==0),
                       colSums(train.y==1), colSums(train.y==2),
                       colSums(train.y==3), colSums(train.y==4),
                       log(colMeans(train.y)), 
                       log(colSums(train.y^2)),
                       log(colSums(train.y^3)),
                       log(colSums(train.y^4)),
                       log(colSums(train.y^5)),
                       log(colSums(train.y^6)),
                       accf.train.y[2:6,], train.y, train.yo, train.y^2, train.yo^2,
                       log(1 + train.y), log(1 + train.yo),
                       train.yad^2, train.yado^2
  )
  }
  
  newlist=list("log_r"=log_r, "log_sig_e"=log_sig_e, "phi"=phi, "sy"=E)
  return(newlist)
}
