GKDR <- function(G, obstemp, theta, focuspara,ss_set) {
  Np = 20000
  Nt = 30000
  N = 2000
  np = 3
  weight_k = 1
  c = 0.75
  
  CV_num = 1000
  # generate training data
   train = simRicker_gkdr(Np,ss_set)
   test = simRicker_gkdr(Nt,ss_set)
  # save(train, test, file='gkdr_test_train.RData')
  #load('gkdr_test_train.RData')
  ##choose set of original summary statistics
  #obstemp=obstemp$E0#1:E0,2:E1,3:E2
  
  test.sy=t(test$sy)
  train.sy = t(train$sy)
  # test.sy.std=apply(test.sy,2,sd)
  train.sy.std=apply(train.sy,2,sd)
  
  train.sy.mean=apply(train.sy,2,mean)
  
  test.sy=sweep(test.sy,2,train.sy.mean,FUN="-")
  test.sy=sweep(test.sy,2,train.sy.std,FUN="/")
  train.sy=sweep(train.sy,2,train.sy.mean,FUN="-")
  train.sy=sweep(train.sy,2,train.sy.std,FUN="/")
  ##
  P = length(obstemp)
  obstemp=(obstemp-train.sy.mean)/train.sy.std
  SumStatDat_v = rep(obstemp, each = Np)
  SumStatDat = matrix(SumStatDat_v, nrow = Np)
  SumStatDatt_v = rep(obstemp, each = Nt)
  SumStatDat_t = matrix(SumStatDatt_v, nrow = Nt)
  # prepare data

  thetas = cbind(train$log_r, train$log_sig_e, train$phi)
  thetast = cbind(test$log_r, test$log_sig_e, test$phi)
  
  sprintf('#sample = %d, dim of X = %d, effective dim = %d\n', N, P, G)
  candx = c(0.25, 0.5, 0.75, 1,2)  # candidates for CV
  candy = c(0.25, 0.5, 1,2)
  eps = 0.00001;
  err_tbl = vector(mode = "numeric",
                   length = (length(candx) * length(candy) * length(eps)))
  err_tbl2 = vector(mode = "numeric",
                   length = (length(candx) * length(candy) * length(eps)))
  err_tbl3 = vector(mode = "numeric",
                   length = (length(candx) * length(candy) * length(eps)))
  #sorting to get bandwidth and training data
  rsse_b = rowMeans((SumStatDat - train.sy) * (SumStatDat - train.sy))
  rsse_t = sort(rsse_b, index.return = TRUE)
  
  Xm = train.sy[rsse_t$ix,]
  Ym = thetas[rsse_t$ix,]
  Y = Ym[1:N,]
  X = Xm[1:N,]
  rtt = rsse_t$x[1:N]
  rsseup = rtt[N]
  W = c * (1 - (rtt / rsseup) ^ 2)
  # Xsd=apply(X,2,sd)
  # Ysd=apply(Y,2,sd)
  # X=sweep(X, 2, obs_sd, FUN="/")
  # Y=sweep(Y, 2, Ysd, FUN="/")
  # X=scale(X, center = FALSE, scale = TRUE)
   Y=scale(Y, center = TRUE, scale = TRUE)
  # train.sy=sweep(train.sy,2,obs_sd,FUN="/")
  # SumStatDat_t=sweep(SumStatDat_t,2,obs_sd,FUN="/")
  sgx0 = MedianDist(X)
  sgy0 = MedianDist(Y)
 # X=X*W
  # parallel validation
  sgx_v = sgx0 * candx
  sgy_v = sgy0 * candy
  tt1 = length(candx) * length(candy) * length(eps)
  tt2 = length(candy) * length(eps)
  for (hhh in 1:tt1) {
    opth = ceiling(hhh / tt2)
    sgx = sgx_v[opth]
    rr = hhh - (opth - 1) * tt2
    optk = ceiling(rr / length(eps))
    sgy = sgy_v[optk]
    opte=rr-(optk-1)*length(eps)

    EPS = eps[opte]
    time.each=proc.time()
    KD = KernelDeriv(X, Y, K, sgx, sgy, EPS, W)
    time.eachend=proc.time()-time.each
    #cat(sprintf("\n time for KernelD in %d round is %f\n",hhh,time.eachend))
    K=8#nScree(x=KD$values)$Components$nkaiser
    B=KD$vectors[,1:K]
    rsse_b = rowMeans((SumStatDat_t%*%B - test.sy%*%B)^2)
    #reesb_rej=rowMeans((SumStatDat_t - test.sy)^2)
    rsse_t = sort(rsse_b, index.return = TRUE)
    #rsse_t_rej = sort(reesb_rej, index.return = TRUE)
    Yt = thetast[rsse_t$ix,]
    #Ytt = thetast[rsse_t_rej$ix,]
    Yt_v = colMeans(Yt[1:100,])
    Yt_v2 = colMeans(Yt[1:10,])
    Yt_v3 = colMeans(Yt[1:5,])
    #Ytt_v = colMeans(Ytt[1:20,])
    err_tbl[hhh] = sum((Yt_v-theta)^2)
    err_tbl2[hhh] = sum((Yt_v2-theta)^2)
    err_tbl3[hhh] = sum((Yt_v3-theta)^2)

    #errr=sum((Ytt_v-theta)^2)
  }
  midx = which.min(err_tbl)
  midx2 = which.min(err_tbl2)
  midx3 = which.min(err_tbl3)
  opth = ceiling(midx3 / tt2)
  sgx = sgx_v[opth]
  rr = midx3 - (opth - 1) * tt2
  optk = ceiling(rr / length(eps))
  sgy = sgy_v[optk]
  opte=rr-(optk-1)*length(eps)

  #set up X,Y and weight again as opimum value
  sgx = sgx0*candx[opth]
  sgy = sgy0*candy[optk]
  EPS = eps[opte]

  KD = KernelDeriv(X, Y, K, sgx, sgy, EPS, W)
  K=8#nScree(x=KD$values)$Components$nkaiser
  para_gkdr = c(sgx, sgy, EPS)
  B=KD$vectors[,1:K]
  
  cat(sprintf("Npop: %d,N: %d", Np, N))
  cat(sprintf(
    "/n K= %d sgx=%f,sgy=%f,eps:%f,Yfocus:%d\n",
    K,
    sgx,
    sgy,
    EPS,
    focuspara
  ))
  return(B)
  }


MedianDist <- function (X) {
  N = length(X[, 1])
  ab = X %*% t(X)
  aa = diag(ab)
  X=matrix(rep(c(aa),times=N),nrow = N)
  Dx=X+t(X)-2*ab
  Dx = Dx - diag(diag(Dx))
  dx = c(Dx)
  dx = dx[dx != 0]
  s = sqrt(median(dx))
}

KernelDeriv <- function(X, Y, K, sgx, sgy, EPS, W){
  dims=dim(X)
  N=dims[1]
  M=dims[2]
  sx2=sgx^2
  sy2=sgy^2
  
  ab=X%*%t(X)
  aa=matrix(rep((diag(ab)),times=N),nrow = N)
  xx=aa+t(aa)-2*ab
  xx[xx<0]=0
  Kx=exp(-xx*(sx2^(-1)))
  
  ab=Y%*%t(Y)
  aa=matrix(rep((diag(ab)),times=N),nrow = N)
  yy=aa+t(aa)-2*ab
  yy[yy<0]=0
  Ky=exp(-yy*(sy2^(-1))) 
  
  F1=solve(Kx+N*EPS*diag(N))%*%Ky%*%solve(Kx+N*EPS*diag(N))
  Mt=array(0,dim=c(M,M,N))
  Xmt=matrix(rep(X,each=N,times=1),nrow = N)
  dim(Xmt)=c(N,N,M)
  Xm=(aperm(Xmt,c(2,1,3))-Xmt)/sx2
  Kxm=rep(Kx,times=M)
  dim(Kxm)=c(N,N,M)
  H=matrix(Xm)*matrix(Kxm)
  dim(H)=c(N,N,M)
  Hm=aperm(H,c(1,3,2))
  for(i in 1:N){
    Mt[,,i]=t(Hm[,,i])%*%F1%*%Hm[,,i]
    Mt[,,i]=matrix(Mt[,,i]*W[i],nrow = M)
  }
  R=rowSums(Mt,dims = 2)
  # Xmt=matrix(rep(X,each=N,times=1),nrow = N)
  # dim(Xmt)=c(N,N,M)
  # Xm=(aperm(Xmt,c(2,1,3))-Xmt)/sx2
  # Kxm=rep(Kx,times=M)
  # H=matrix(Xm,nrow=N)*matrix(Kxm,nrow=N)
  # Hm=t(H)%*%H
  # #F1t=matrix(rep(t(F1),times=M),ncol = N,byrow = TRUE)
  # #F1tt=t(matrix(rep(F1t),times=M,ncol = N,byrow = TRUE))
  # dim(Hm)=c(N,M,N,M)
  # Hm=aperm(Hm,c(1,3,2,4))
  # dim(Hm)=c(N*N,M,M)
  # F1=rep(F1,times=M*M)
  # dim(F1)=c(N*N,M,M)
  # R=colSums(Hm*F1,dims=1)
  # dim(R)=c(M,M)
  
  V=eigen(R, symmetric = TRUE)
  
}
