function [B] = GKDR_fix(K,c)
%GKDR Summary of this function goes here
%   Detailed explanation goes here
% Bm=fopen('Bmatrix_s.txt','a+');
N = para('training_number');
np = para('np');
M = para('n_summarystatics');

X = zeros(N,M);   %t:N,ns:M
Y_t=zeros(N,np);
W=zeros(N,1);
load('tY.mat');
load('train.mat');
local_band=para('local_b');
sg_f=[26.62,7];

K=2;

%%%%%%%%%%%%%%%%%%%%%%%%%
%Localized
SumStatDat = [28 5 2 3 0 0 5];

%%%%%%%%%%%%%%%%%%%%%%%%%
%Localized
ss=0;
% iii=1;
% while ss<200
%     iii=iii+1;
for iii=1:N
    SumStatSim_f=msdata(iii,:);
    SumStatDat_f=SumStatDat;
    n1 = norm(minus(SumStatSim_f,SumStatDat_f));
    if(SumStatSim_f==0)
        continue;
    end
    if(n1 <= local_band)
        ss = ss + 1;
        X(ss,:) = msdata(iii,:);
        Y_t(ss) = Y(iii);
       %W(ss) = c*(1-(n1/local_band)^2);
       W(ss)=1-n1/local_band;
    end
end

fprintf('training:%3d\n',ss);
%%%%%%%%%%%%%%%
%Update our buffer
N=ss;
X = X(any(X,2),:);   %t:N,ns:M
Y_t=Y_t(any(Y_t,2),:);

W=W(W~=0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Prepare Parameters
Y = Y_t;
sgx=sg_f(1);
sgy=sg_f(2);
EPS=0.001;
[B]=KernelDeriv(X,Y,K,sgx,sgy,EPS,W);


end


