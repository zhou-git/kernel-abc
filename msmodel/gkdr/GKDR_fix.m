function [B] = GKDR_fix(K,c,obsdata)
%GKDR Summary of this function goes here
%   Detailed explanation goes here
% Bm=fopen('Bmatrix_s.txt','a+');
N = para('training_number');
np = para('np');
M = para('n_summarystatics');

load('tdata.mat');
load('tY.mat');
load('ms_20_1M_mad.mat');

% X = zeros(N,M);   %t:N,ns:M
% Y_t=zeros(N,np);

tnum=30000;%para('tdata_num');

A=diag(1./ms_20_1M_mad);
SumStatDat = obsdata;

% %%%%%%%%%%%%%%%%%%%%%%%%%
% ss=0;
% iii=0;
% while ss<N
%     iii=iii+1;
%     SumStatSim_f=tdata(iii,:);
%     SumStatDat_f=SumStatDat;
%     n1 = minus(SumStatSim_f,SumStatDat_f)*A*minus(SumStatSim_f,SumStatDat_f)';%norm(minus(SumStatSim_f,SumStatDat_f));
%     if(SumStatSim_f==0)
%         continue;
%     end
%     if(n1 < local_band)
%         ss = ss + 1;
%         X(ss,:) = tdata(iii,:);
%         Y_t(ss) = Y(iii);
%        % W(ss) = c*(1-(n1/local_band)^2);
%        W(ss)=1-n1/local_band;
%        if W(ss)==0
%            sds=0;
%        end
%     end
% end
rsse_b=zeros(1,N);
rsse_t=zeros(1,N);


for i=1:tnum;%(round(N/size_data))

    SumStatSim = tdata(i,:);%[SumStatMaxSim SumStatMinSim SumStatQuntSim]';
    rsse_b(i) = minus(SumStatSim,SumStatDat)*A*minus(SumStatSim,SumStatDat)'; 
end
% fprintf('X %f %f %f %f %f %f %f\n Y %f\n',sum(X)./ss,sum(Y_t)./ss);
% fprintf('training:%3d total %d\n',ss,iii);
%%%%%%%%%%%%%%%
%Update our buffer
%N=ss;
[rsse_t,I]=sort(rsse_b);
Ym=tY(I);
Xm=tdata(I,:);
W=1.-rsse_t(1:N)./rsse_t(N);
W(W==0)=0.00001;
W2 = c.*(1.-(rsse_t(1:N)./rsse_t(N+1)).^2);
Y_t=Ym(1:N);
X=Xm(1:N,:);


% fprintf('training:%3d\n',N);
%%%%%%%%%%%%%%%
%Update our buffer
% N=ss;
% X = X(any(X,2),:);   %t:N,ns:M
% Y_t=Y_t(any(Y_t,2),:);

W=W2';%ones(N,1);%W2';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Prepare Parameters
Y = Y_t;
sgx=17;%sg_f(1);
sgy=2;%sg_f(2);
EPS=0.0001;
[B]=KernelDeriv(X,Y,K,sgx,sgy,EPS,W);
fprintf('B finished\n');
end


