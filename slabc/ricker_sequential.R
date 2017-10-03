#!/usr/bin/env Rscript
#install.packages("slabc_0.0-5.tar.gz",repos=NULL)
rm(list=ls())
library(slabc)
library(nFactors)
setMKLthreads(32)
sink("./output/output.txt")
#load observation
load("rickerFIX.Rdata")
#warning: places to change the summary statistics set including : ricker_se*.R, GKDR.R(two function),and internal.R
#
ricker_prior=list(c("unif",0,10),c("unif",log(0.1),0),c("unif",0,30)) # a uniform prior distribution between 0 and 1
## define the targeted summary statistics
theta=c(obs.log_r[1],obs.log_sig_e[1],obs.phi[1])
##################################################
#Parameters
ss_set=1
alpha_delmo=0.5
tolerance=0.01
g_flag=1
nbsimu=5000
#############################################
sum_stat_obs=list("E0"=obs.sy.E0[,1],"E1"=obs.sy.E1[,1],"E2"=obs.sy.E2[,1])
if(ss_set==1){
  sum_stat_obs = sum_stat_obs$E1
  sd_obs=sd_E1
}else if(ss_set==0){
  sum_stat_obs = sum_stat_obs$E0
  sd_obs=sd_E0
}else{
  sum_stat_obs = sum_stat_obs$E2
  sd_obs=sd_E2
}


##
time.begin=proc.time()

cat("to perform the Del Moral et al. (2012) s method: ricker model\n")
print(paste("ss_set=",ss_set,"tolerance=",tolerance,"g_flag=",g_flag,"nb_simul=",nbsimu,sep = " "))
if(g_flag==1){
  
  G=8
  focuspara=0
  pB=GKDR(G, sum_stat_obs, theta, focuspara,ss_set)
  save(pB,file = "B_matrixE1.RData")
#  load("Bm.RData")
  time.GKDR.full=proc.time()-time.begin
  time.smc.begin=proc.time()
  cat(sprintf("\n time for GKDR full is %f\n", sum(time.GKDR.full)))
  ABC_Delmoral<-ABC_sequential(method="Delmoral", model=simuRicker, prior=ricker_prior,
                               nb_simul=nbsimu, summary_stat_target=sum_stat_obs, 
                               tolerance_target=tolerance,
                               g_flag,pB,
                               verbose = FALSE,
                               use_seed = TRUE
                               #n_cluster = 32
  )
  time.SMC.full= proc.time()-time.smc.begin
  cat(sprintf("\n time for Sequential MC is %f", sum(time.SMC.full)))
  save(ABC_Delmoral,file = "result_gkdr.RData")
 
}else if(g_flag==0){
  cat(sprintf("SMC for ricker\n"))
  time.SMC.begin=proc.time()
  ABC_Delmoral<-ABC_sequential(method="Delmoral", model=simuRicker, prior=ricker_prior,
                               nb_simul=nbsimu, summary_stat_target=sum_stat_obs, tolerance_target=tolerance,g_flag,pB=NULL,verbose = FALSE,
                               
                               use_seed = TRUE)
  time.SMC.full = proc.time()-time.SMC.begin
  cat(sprintf("\n time for Sequential MC is %f", sum(time.SMC.full)))
  save(ABC_Delmoral,file = "result_smc.RData")
  
}else{
  cat(sprintf("Semi for Ricker\n"))
  time.semi.begin=proc.time()
  pB=SEMI(sum_stat_obs, theta,ss_set,sd_obs)
  time.semi.end=proc.time()-time.semi.begin
  time.smc.begin=proc.time()
  cat(sprintf("time for semi is %f \n",sum(time.semi.end)))
  ABC_Delmoral<-ABC_sequential(method="Delmoral", model=simuRicker, prior=ricker_prior,
                               nb_simul=nbsimu, summary_stat_target=sum_stat_obs, 
                               tolerance_target=tolerance,
                               g_flag,pB$B,
                               verbose = FALSE,
                               #n_cluster = 32,
                               use_seed = TRUE) 
  time.semi.full= proc.time()-time.smc.begin
  cat(sprintf("\n time for Sequential MC is %f", sum(time.semi.full)))
  save(ABC_Delmoral,file = "result_semi.RData")
}

sink()
ABC_Delmoral
