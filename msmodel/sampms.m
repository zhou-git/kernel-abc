function [q_out,mser_out,mean_out] = sampms(theta,c,delta)
%Rejection method
%
Acc = 1000;
N=16000;
M=7;

LLR = para('LLR');
G=para('effective_dimension');
diff_SS_buffer = zeros(M,Acc);    %%Memory run out, do not need so large matrix for weights
weight_buffer = zeros(Acc,1);
theta_buffer = zeros(1,Acc);

%%%%%%%%%%%%%%%%%%%%%%%%%
%Prepare Data

SumStatDat = [28 5 2 3 0 0 5];
%SumStatDat = grtdata;
starttime = cputime;
%%%%%%%%%%%%%%%%%%%%%%%%%%
%Simulation
load('Y.mat');%sdata 1xe7
load('msdata.mat');%stheta
i=0;
j=0;
while i<Acc && j<N
    j=j+1;
    SumStatSim = msdata(j,:);%[SumStatMaxSim SumStatMinSim SumStatQuntSim]';
    n1 = norm(minus(SumStatSim,SumStatDat));
    if(n1 <= delta)
        i = i + 1;
        theta_buffer(i) = Y(j);
        diff_SS_buffer(:,i) = SumStatSim-SumStatDat;
        weight_buffer(i,1) = c/delta*(1-(n1/delta)^2);
    end
   % simss(:,j) = SumStatSim-SumStatDat;
end
    
fprintf('number of accepted samples=%3d\n',i);

if(LLR)
    mean_out = LocalLinearRgs(theta_buffer, diff_SS_buffer,weight_buffer,i,M);
    mser_out = ((mean_out - theta).^2)./(theta.^2);
    q_out = i/j;
else
    thetam = sum(theta_buffer);
    mean_out = thetam'/i;              %%use mean as point estimation
    mser_out = ((mean_out - theta).^2)./(theta.^2);
    q_out = i/j;
end  

endtime = cputime - starttime;
fprintf('"Simulation(rejection), datasize = %3.0f"\n',Acc);
fprintf('Time = %3.4f\n', endtime);

end
