library(R.matlab)
writeMat(con="Y.mat", log_r=log_r, log_sig_e = log_sig_e, phi=phi, sYE0=sY.E0, sYE1=sY.E1, Y=Y, matVersion='5')
sYE2_1=sY.E2[1:150,]
sYE2_2=sY.E2[151:300,]
sYE2_3=sY.E2[301:424,]
writeMat(con="sYE2_1.mat",sYE2_1=sYE2_1,matVersion='5')
writeMat(con="sYE2_2.mat",sYE2_2=sYE2_2,matVersion='5')
writeMat(con="sYE2_3.mat",sYE2_3=sYE2_3,matVersion='5')
writeMat(con="train.mat",tlog_r=train.log_r, tlog_sig_e = train.log_sig_e, tphi=train.phi, tsyE0=train.sy.E0, tsyE1=train.sy.E1, tsyE2=train.sy.E2, ty=train.y, matVersion='5')
writeMat(con="yobs3.mat", obs_log_r=obs.log_r, obs_log_sig_e=obs.log_sig_e, obs_phi=obs.phi, obs_syE0=obs.sy.E0, obs_syE1=obs.sy.E1, obs_syE2=obs.sy.E2,obs_y=obs.y, matVersion='5')
