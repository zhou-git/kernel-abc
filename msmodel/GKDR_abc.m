function [B_out] = GKDR_abc(G,c)
%GKDR Summary of this function goes here
%   Detailed explanation goes here
% Bm=fopen('Bmatrix_s.txt','a+');
N = para('training_number');
np = para('np');
M = para('n_summarystatics');

K = G;
X = zeros(N,M);   %t:N,ns:M
Y_t=zeros(N,np);
W=zeros(N,1);
X_b=zeros(N,M);
Y_b=zeros(N,np);

local_band=para('local_b');
load('tY.mat');
load('train.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%
%Simulation     1,no rejection 2.with rejection
%N: number of training data for all parameters(different B)
%1:
starttime=cputime;

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
       % W(ss) = c*(1-(n1/local_band)^2);
       W(ss)=1-n1/local_band;
       if W(ss)==0
           sds=0;
       end
    end
end

%%%%%%%%%%%%%%%
%Update our buffer
N=ss;
X = X(any(X,2),:);   %t:N,ns:M
Y_t=Y_t(any(Y_t,2),:);
W=W(W~=0);

fprintf('Cross validation ...\n');
fprintf('#sample = %d, dim of X = %d, effective dim = %d\n', N, M, K);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%Prepare Parameters

    Y = Y_t;
   % Y=Y_all(:,1);
   % K =1;
    sgx0=MedianDist(X);   % Basic value for bandwidth for X
    sgy0=MedianDist(Y);   % Basic value for bandwidth for Y
    NCV=5;      % Number of cross-validation
    DEC=10;     % Number of decreasing dimensions for gKDR-i;

    B0=zeros(M,1);
    B0(1,1)=1;      % B: Projection matrix onto the effective direction

    candx=[0.25 0.5 0.75 1 2];  % candidates for CV
    candy = [0.25 0.5 0.75 1 2];
    eps = [0.000001];%[0.00001];

    % For cross-validation
    ridx=randperm(N);  % random order 
    Xr=X(ridx,:);
    Yr=Y(ridx,:);  
    Wr=W(ridx,:);
    lx=ceil(N/NCV);
    ei=cumsum(lx.*ones(1,NCV),2);  
    si=ei-(lx-1).*ones(1,NCV);
    ei(NCV)=N;       % si: staring idx, ei: ending idx
    err_tbl=zeros(length(candx)*length(candy)*length(eps), NCV);

    for h=1:length(candx)
        sgx=sgx0*candx(h);% five different values for CV
        for k=1:length(candy)
          sgy=sgy0*candy(k);            
          for ll=1:length(eps)
            EPS = eps(ll);
            for i=1:NCV
                ri=si(i):ei(i);
                Xe=Xr; Ye=Yr; We=Wr;
                Xe(ri,:)=[];
                Ye(ri,:)=[];    % Xe, Ye: trainig sample for CV
                We(ri,:)=[];
                Xt=Xr(ri,:);
                Yt=Yr(ri,:);    % Xt, Yt: test sample for CV
                [B, t]=KernelDeriv(Xe,Ye,K,sgx,sgy,EPS,We);

            % kNN regression for CV
                nnidx=knnsearch(Xe*B,Xt*B, 'K', 5, 'NSMethod', 'kdtree');

                Yo=zeros(length(ri),length(Y(1,:)));
                for j=1:length(ri)
                    ii=nnidx(j,:);
                    Yo(j,:)=sum(Ye(ii',:),1)./5;
                end

                dd=Yt-Yo;      
                err_tbl((h-1)*length(candy)*length(eps)+(k-1)*length(eps)+ll,i)=sum(sum(dd.*dd,1),2)./length(ri);  
                %5 means divide the data to 5 parts. 
            end
          end
        end
    end
    
%     fprintf('Validation Done\n');
    [c, midx]=min(mean(err_tbl,2));
    opth=ceil(midx/(length(candy)*length(eps)));
    rr=midx-(opth-1)*length(candy)*length(eps);
    optk=ceil(rr/length(eps));
    opte=mod(rr,length(eps));
    if opte==0
        opte=length(eps);
    end

    % Parameter
    sgx=sgx0*candx(opth);
    sgy=sgy0*candy(optk);
    EPS=eps(opte);
    [B, t]=KernelDeriv(X,Y,K,sgx,sgy,EPS,W);

%   err=sqrt(trace(B0*B0'*(eye(M)-B*B'))/trace(B0'*B0));

    B_out=B;
%   fprintf('NCV=%d, DEC=%d, K=%d\n', NCV, DEC);
    fprintf('sgx = %.3f  sgy = %.3f eps=%e\n', sgx,sgy,EPS);
%     fprintf('Extimation Error: %f\n', mean(err));

% fprintf(Bm,'\nB_out matrix\n%3f\n',B);
% fclose(Bm);
% disp(B_out);
endtime = cputime - starttime;
fprintf('Time = %3.4f\n', endtime);

end

