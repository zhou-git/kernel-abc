function [q,mser,thetam] = RejGKDR(data,theta,c,delta,B)
%Rejection method
%
Acc = 5000;
N=30000000;
NQ=para('num_quantile');
LLR = para('LLR');
G=para('effective_dimension');

%%Buffer
diff_SS_buffer = zeros(G,N);    %%Memory run out, do not need so large matrix for weights
weight_buffer = zeros(N,1);
theta_buffer = zeros(N,3);

%%Prepare Data
SumStatMaxDat = max(data);  %summary statistics 1:max
SumStatMinDat = min(data);  %...                2:min
SumStatQuntDat = quantile(data,NQ); %           3:quantile of 3
SumStatDat = [SumStatMaxDat SumStatMinDat SumStatQuntDat];

starttime = cputime;

fprintf('Bandwith %3.3f    ',delta);
%Calculating
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GKDR
i=0;
j=0;
while i<Acc && j<N
    j=j+1;
    Pri1 = 10*rand(1,1); 
    Pri2 = Pri1 + 10*rand(1,1);
    Pri3 = 1/3*rand(1,1);
    YSim = simulationMG1([Pri1 Pri2 Pri3]);
    SumStatQuntSim = quantile(YSim,NQ);
    SumStatMaxSim = max(YSim);
    SumStatMinSim = min(YSim);
    SumStatSim = [SumStatMaxSim SumStatMinSim SumStatQuntSim];
    SumStatSim_f=SumStatSim*B;
    SumStatDat_f=SumStatDat*B;
    n1 = norm(minus(SumStatSim_f,SumStatDat_f));
    if(n1 <= delta)
        i = i + 1;
        theta_buffer(i,:) = [Pri1 Pri2 Pri3];
        diff_SS_buffer(:,i) = SumStatSim_f-SumStatDat_f;
        weight_buffer(i,1) = c/delta*(1-(n1/delta)^2);
    end
end
if(LLR)
    thetam = LocalLinearRgs(theta_buffer, diff_SS_buffer,weight_buffer,i,G);
    mser = ((thetam - theta).^2)./(theta.^2); 
    q = i/j;
else
    thetam = sum(theta_buffer,1);
    thetam = thetam/i;              %%use mean as point estimation
    mser = ((thetam' - theta).^2)./(theta.^2);
    q = i/j;
end 

endtime = cputime - starttime;
fprintf('Time = %3.4f\n', endtime);

end
