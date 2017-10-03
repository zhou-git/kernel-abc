SEMI = function(sum_stat_obs,theta,ss_set,sd_obs){
  Np = 20000
  N = 3000
  np = 3
  
  # generate training data
  train = simRicker_gkdr(Np,ss_set)
  train.sy=t(train$sy)
  train.sy_mean=apply(train.sy,2,mean)
  train.sy_sd=apply(train.sy,2,sd)
  train.sy=sweep(train.sy,2,train.sy_mean,FUN="-")
  train.sy=sweep(train.sy,2,train.sy_sd,FUN="/")
  ##
  obstemp=(sum_stat_obs-train.sy_mean)/train.sy_sd
  P = length(obstemp)
  SumStatDat_v = rep(obstemp, each = Np)
  SumStatDat = matrix(SumStatDat_v, nrow = Np)
  # prepare data
  
  thetas = cbind(train$log_r, train$log_sig_e, train$phi)

  rsse_b = rowMeans((SumStatDat - train.sy) * (SumStatDat - train.sy))
  rsse_t = sort(rsse_b, index.return = TRUE)
  
  Xm = train.sy[rsse_t$ix,]
  Ym = thetas[rsse_t$ix,]
  Y = Ym[1:N,]
  X = Xm[1:N,]
  rtt = rsse_t$x[1:N]
  rsseup = rtt[N]
  
  
  
  Y=scale(Y, center = TRUE, scale = TRUE)
  EPS=0.0001
  # Z=as.data.frame(cbind(Y,X))
  # reg=lm(Y~.,Z)
  X=cbind(1,X)
  d=length(theta)
  
  B0 <- c()
  B <- matrix(nrow=P, ncol=d)
  for (i in 1:d) { 
    reg <- lm.fit(X,Y[,i])
    B0[i] <- reg$coefficients[1]
    B[,i] <- reg$coefficients[-1]
  }
  if(any(is.na(B))) warning("Linear regression ill-conditioning")
  B[is.na(B)]=0
  return(list(B0=B0,B=B))
}