library(slabc)
# library(parallel)
#library(R.matlab)
set.seed(2);
# cl <- makeForkCluster(16 )
snum<-1e+02
##Parameter setting 
n.total <- snum*100
burn.in <- 50
n.t     <- 50

obs.log_r     <- c(rep(3.8,n.total))#runif(n.total, min = 0, max = 20)#
obs.log_sig_e <- runif(n.total, min = log(0.1), max = 0)
obs.sig_e     <- exp(obs.log_sig_e)
obs.phi       <- c(rep(10,n.total))#runif(n.total, min = 0, max <- 30)

e <- matrix(rnorm((burn.in+n.t) * n.total), burn.in + n.t, n.total)
u <- runif(n.t * n.total) 
n.total <- snum
## simulate data ....
obs.yt <- ricker_abc(obs.phi, obs.log_r, obs.sig_e, e, u, burn.in)

obs.yt[obs.yt==0]<-NA

ttable <-colSums(!is.na(obs.yt))>20
obs.ytm <- obs.yt[,ttable]

obs.log_rf <- obs.log_r[ttable]
obs.log_sig_ef<- obs.log_sig_e[ttable]
obs.sig_ef<- obs.sig_e[ttable]
obs.phif<- obs.phi[ttable]

obs.ytm[is.na(obs.ytm)] <- 0
obs.y<-obs.ytm[,1:n.total]

obs.log_r <- obs.log_rf[1:n.total]
obs.log_sig_e<- obs.log_sig_ef[1:n.total]
obs.sig_e<- obs.sig_ef[1:n.total]
obs.phi<- obs.phif[1:n.total]
## Now assemble the relevant statistics
acf.obs.y  <- sl.acf(obs.y,max.lag=5,0)
accf.obs.y <- sl.acf(obs.y,max.lag=5,1)
b0.obs.y <- nlar(obs.y^.3,lag=c(1,1),power=c(1,2))
obs.yad <- obs.y[2:n.t, ] - obs.y[1:n.t-1, ]
b1.obs.y <- order.dist(obs.y,obs.yad,np=3,diff=1)  

obs.yo <- apply(obs.y, 2, sort)
obs.yad <- obs.y[2:n.t, ] - obs.y[1:n.t-1, ]  #time difference
obs.yado <- apply(obs.yad, 2, sort)

## combine the statistics...
obs.sy.E0 <- rbind(acf.obs.y,
               b0.obs.y, b1.obs.y,
               colMeans(obs.y),
               colSums(obs.y==0)
)

obs.sy.E1 <- rbind(acf.obs.y, b0.obs.y, b1.obs.y,
               colMeans(obs.y), colSums(obs.y==0),
               colSums(obs.y==1), colSums(obs.y==2),
               colSums(obs.y==3), colSums(obs.y==4),
               log(colMeans(obs.y)), 
               log(colSums(obs.y^2)),
               log(colSums(obs.y^3)),
               log(colSums(obs.y^4)),
               log(colSums(obs.y^5)),
               log(colSums(obs.y^6)),
               accf.obs.y[2:6,]
)
obs.sy.E2 <- rbind(acf.obs.y, b0.obs.y, b1.obs.y,
               colMeans(obs.y), colSums(obs.y==0),
               colSums(obs.y==1), colSums(obs.y==2),
               colSums(obs.y==3), colSums(obs.y==4),
               log(colMeans(obs.y)), 
               log(colSums(obs.y^2)),
               log(colSums(obs.y^3)),
               log(colSums(obs.y^4)),
               log(colSums(obs.y^5)),
               log(colSums(obs.y^6)),
               accf.obs.y[2:6,], obs.y, obs.yo, obs.y^2, obs.yo^2,
               log(1 + obs.y), log(1 + obs.yo),
               obs.yad^2, obs.yado^2
)
sd_E0=apply(obs.sy.E0,1,sd)
mean_E1=apply(obs.sy.E1,1,mean)
sd_E1=apply(obs.sy.E1,1,sd)
mean_E0=apply(obs.sy.E0,1,mean)
sd_E2=apply(obs.sy.E2,1,sd)
mean_E2=apply(obs.sy.E2,1,mean)
#writeMat(con="yobsFix.mat", obs_log_r=obs.log_r, obs_log_sig_e=obs.log_sig_e, obs_phi=obs.phi, obs_syE0=obs.sy.E0, obs_syE1=obs.sy.E1, obs_syE2=obs.sy.E2,obs_y=obs.y, matVersion='5')
save(obs.log_r, obs.log_sig_e, obs.phi, obs.y, obs.sy.E0, obs.sy.E1, obs.sy.E2,sd_E2,sd_E1,sd_E0,mean_E1,mean_E2,mean_E0, file = "rickerFIX.Rdata")
#save(obs.log_r, obs.log_sig_e, obs.phi, obs.y, obs.sy.E0, obs.sy.E1, obs.sy.E2, file = "rickerRND.Rdata")
