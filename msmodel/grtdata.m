function [Sumstat_buffer]=grtdata
N=16000;
M=7;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%Simulation
load('tY.mat');%sdata 1xe7
load('tmsdata.mat');%stheta
load('data.mat');
delta=20;
for j=1:N
    SumStatSim = msdata(j,:);%[SumStatMaxSim SumStatMinSim SumStatQuntSim]';
    n1 = norm(minus(SumStatSim,data));
    if(n1 <= delta)
        theta_buffer = Y(j);
        Sumstat_buffer = SumStatSim;
        delta=n1;
    end
   % simss(:,j) = SumStatSim-SumStatDat;
end

save theta theta_buffer; 
end