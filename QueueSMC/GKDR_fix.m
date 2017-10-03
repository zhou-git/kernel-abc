function [B] = GKDR_fix(data,K,c,index)
%GKDR Summary of this function goes here
%   Detailed explanation goes here
% Bm=fopen('Bmatrix_s.txt','a+');
np = para('np');
M = para('n_summarystatics');
NQ= para('num_quantile');
Acc=para('Accepted_trainingsample');
N=10000000;
X = zeros(N,M);   %t:N,ns:M
Y_t=zeros(N,np);
W=zeros(N,1);

local_band=para('local_b');
sg_f=[6.333 1.6;9.221 2.402;3.7 0.04];

%%%%%%%%%%%%%%%%%%%%%%%%%
%Prepare Data
SumStatMaxDat = max(data);  %summary statistics 1:max
SumStatMinDat = min(data);  %...                2:min
SumStatQuntDat = quantile(data,NQ); %           3:quantile of 3
SumStatDat = [SumStatMaxDat SumStatMinDat SumStatQuntDat];

%%%%%%%%%%%%%%%%%%%%%%%%%%
%Simulation     1,no rejection 2.with rejection
%N: number of training data for all parameters(different B)
%1:
starttime=cputime;
ss=0;
j=0;
while ss<Acc
    j=j+1;
    Pri1 = 10*rand(1,1); 
    Pri2 = Pri1 + 10*rand(1,1);
    Pri3 = 1/3*rand(1,1);
    YSim = simulationMG1([Pri1 Pri2 Pri3]);
    SumStatQuntSim = quantile(YSim,NQ);
    SumStatMaxSim = max(YSim);
    SumStatMinSim = min(YSim);
    SumStatSim = [SumStatMaxSim SumStatMinSim SumStatQuntSim];
    n1 = norm(minus(SumStatSim,SumStatDat));
    if(n1 <= local_band)
        ss = ss + 1;
        X(ss,:) = SumStatSim;
        Y_t(ss,:) = [Pri1 Pri2 Pri3];
        W(ss)=1-n1/local_band;
        % W(ss) = c*(1-(n1/local_band)^2);
    end
end
fprintf('training:%3d    ',ss);
%%%%%%%%%%%%%%%
%Update our buffer
%N=ss;
X = X(any(X,2),:);   %t:N,ns:M
Y_t=Y_t(any(Y_t,2),:);
W=W(W~=0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Prepare Parameters
Y = Y_t(:,index);


sgx=sg_f(index,1);
sgy=sg_f(index,2);
EPS=0.000001;
[B]=KernelDeriv(X,Y,K,sgx,sgy,EPS,W);

endtime = cputime - starttime;
fprintf('Time = %3.4f\n', endtime);

end


