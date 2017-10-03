function [q,mser,thetam] = RejGKDR(theta,c,delta,B)
%Rejection method
%
Acc = 3000;
N=1000000;
LLR = para('LLR');
G=para('effective_dimension');
G=2;
%%Buffer
diff_SS_buffer = zeros(G,N);    %%Memory run out, do not need so large matrix for weights
weight_buffer = zeros(N,1);
theta_buffer = zeros(N,1);


SumStatDat = [28 5 2 3 0 0 5];

%Load simulated data
load('msdata.mat');%sdata 1xe7
load('Y.mat');%stheta
starttime = cputime;

fprintf('Bandwith is %3.3f\n',delta);
%Calculating
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GKDR
i=0;
j=0;
while i<Acc && j<N
    j=j+1;
    SumStatSim_f=msdata(j,:)*B;
    SumStatDat_f=SumStatDat*B;
    n1 = norm(minus(SumStatSim_f,SumStatDat_f));
    if(n1 <= delta)
        i = i + 1;
        theta_buffer(i,:) = Y(j);%[Pri1; Pri2; Pri3];
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
fprintf('Accpt=%d\n',i);

end
