function [B_out] = GKDR_single(K,fix,SSdata,i,obstheta)
%GKDR Summary of this function goes here
%   Detailed explanation goes here

Np =10000;
N=3000;

c=0.75;
load('obsdata.mat');
load('train.mat');

% X = zeros(N,M);   %t:N,ns:M
% Y_t=zeros(N,np);
tsyE=tSS; 
tSS_gkdr = tSS;
%tsyE=sort(tY,2)';
thetas=tPri;
[tnum,M]=size(tSS_gkdr);

A=diag(1./var(tsyE,0,1));%% distance is squared

%%%%%%%%%%%%%%%%%%%%%%%%%%
%Simulation     1,no rejection 2.with rejection
%N: number of training data for all parameters(different B)
%1:
starttime=cputime;
SumStatDat=SSdata;

% tnum=60000;
rsse_b=zeros(1,tnum);


for in=1:tnum;%(round(N/size_data))

    SumStatSim = tsyE(in,:);%[SumStatMaxSim SumStatMinSim SumStatQuntSim]';
    rsse_b(in) = minus(SumStatSim,SumStatDat)*A*minus(SumStatSim,SumStatDat)'; 
end
[rsse_t,I]=sort(rsse_b);
Ym=thetas(I,:);
%Xm=tsyE(:,I)';
Xm=tSS_gkdr(I,:);

Y1=Ym(1:Np,:);
X1=Xm(1:Np,:);

% 
%  Nn=randperm(Np);
%  X2=X1(Nn,:);
%  Y2=Y1(Nn,:);
%  Ytt=Y2(1:N,:);
%  Xtt=X2(1:N,:);



% Xttm=mean(X2);
Xstd=std(X1);
Xstd(Xstd==0)=1;
Ystd=std(Y1);
X=X1;%./repmat(Xstd,N,1);
%Yttm=mean(Y2);
Y_t=Y1;%./repmat(Ystd,N,1);


focuspara=0;

if(focuspara==0)
    Y=Y_t;
    obsY=obstheta';
elseif(focuspara==1)
    Y=Y_t(:,1);
    obsY=obstheta(1);
elseif(focuspara==2)
    Y=Y_t(:,2);
    obsY=obstheta(2);
elseif(focuspara==3)
    Y=Y_t(:,3);
    obsY=obstheta(3);
end
%fprintf('training:%3d\n',N);
%%%%%%%%%%%%%%%
%Update our buffer
% N=ss;
% X = X(any(X,2),:);   %t:N,ns:M
% Y_t=Y_t(any(Y_t,2),:);
if(fix==0)
    fprintf('Cross validation ...\n');
    fprintf('#sample = %d, dim of X = %d, effective dim = %d\n', N, M,K);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%Prepare Parameters
    %Y = Y_t;
    % Y=Y_all(:,1);
    % K =1;
    sgx0=MedianDist(X);   % Basic value for bandwidth for X
    sgy0=MedianDist(Y);   % Basic value for bandwidth for Y
    NCV=5;      % Number of cross-validation
    DEC=10;     % Number of decreasing dimensions for gKDR-i;

%     B0=zeros(M,1);
%     B0(1,1)=1;      % B: Projection matrix onto the effective direction

    candx=[0.25 0.5 0.75 1 2];  % candidates for CV
    candy = [0.25 0.5 0.75 1 2];
    eps = [0.0001 0.001 0.00001];%[0.00001];

    % For cross-validation
    ridx=randperm(N);  % random order 
    Xr=X(ridx,:);
    Yr=Y(ridx,:);  
    
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
                Xe=Xr; Ye=Yr; 
                Xe(ri,:)=[];
                Ye(ri,:)=[];    % Xe, Ye: trainig sample for CV

                Xt=Xr(ri,:);
                Yt=Yr(ri,:);    % Xt, Yt: test sample for CV
                %[B]=KernelDeriv_for(Xe,Ye,K,sgx,sgy,EPS,We);
                %for loop kdr
                NCVi=5;
                [Nt,~]=size(Xe);
                lxt=ceil(Nt/NCVi);
                eit=cumsum(lxt.*ones(1,NCVi),2);  
                sit=eit-(lxt-1).*ones(1,NCVi);
                eit(NCVi)=Nt;       % si: staring idx, ei: ending idx

                BM=zeros(M,M,3);
                if(NCVi>1)
                    for ii=1:NCVi
                        Xtem=Xe(sit(ii):eit(ii),:);
                        Ytem=Ye(sit(ii):eit(ii),:);
                     
                        [~,BM(:,:,ii)]=KernelDeriv_single(Xtem,Ytem,K,sgx,sgy,EPS,SSdata,obsY);%,NDIV);
                    end
                    %fprintf('done!\n');
                    BMtotal=sum(BM,3);
                    [B,~]=eigs(BMtotal,K);
                else
                    [B]=KernelDeriv_single(X,Y,K,sgx,sgy,EPS,SSdata,obsY);
                end
                

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
    %[B,~]=KernelDeriv_for(X,Y,K,sgx,sgy,EPS,W);
    NCVf=5;
    [Nf,~]=size(X);
    lxf=ceil(Nf/NCVf);
    eif=cumsum(lxf.*ones(1,NCVf),2);  
    sif=eif-(lxf-1).*ones(1,NCVf);
    eif(NCVf)=Nf;       % si: staring idx, ei: ending idx

    BF=zeros(M,M,3);
    if(NCVi>1)
        for ii=1:NCVi
            Xtem2=X(sif(ii):eif(ii),:);
            Ytem2=Y(sif(ii):eif(ii),:);
           
          [~,BM(:,:,ii)]=KernelDeriv_single(Xtem2,Ytem2,K,sgx,sgy,EPS,SSdata,obsY);%,NDIV);
        end
        BFtotal=sum(BF,3);
        [B,~]=eigs(BFtotal,K);
    else    
        [B]=KernelDeriv_single(X,Y,K,sgx,sgy,EPS,SSdata,obsY);%,NDIV);

    end
%     err=sqrt(trace(B0*B0'*(eye(M)-B*B'))/trace(B0'*B0));

    B_out=B;
    %   fprintf('NCV=%d, DEC=%d, K=%d\n', NCV, DEC);
    fprintf(' Np: %d,N: %d',Np,N);
    fprintf(',G= %d sgx=%f,sgy=%f,eps:%f,Yfocus:%d\n',K,sgx,sgy,EPS,focuspara);

%     fprintf('Extimation Error: %f\n', mean(err));
else
%     Y = Y_t(:,1);
    NCVi=1;
    sgx=25.8;    %  E3:10,2  E2   E1 13.5 0.9/0.5
    sgy=8.0;%2;
    EPS=0.001;
    lx=ceil(N/NCVi);
    ei=cumsum(lx.*ones(1,NCVi),2);  
    si=ei-(lx-1).*ones(1,NCVi);
    ei(NCVi)=N;       % si: staring idx, ei: ending idx
    if(i==1)
        fprintf(' Np: %d,N: %d',Np,N);
        fprintf(',G= %d sgx=%f,sgy=%f,eps:%f,Yfocus:%d\n',K,sgx,sgy,EPS,focuspara);
    end
    BM=zeros(M,M,3);
    if(NCVi>1)
        for ii=1:NCVi
            Xtem=X(si(ii):ei(ii),:);
            Ytem=Y(si(ii):ei(ii),:);
            Wt=W(si(ii):ei(ii));
            [~,BM(:,:,ii)]=KernelDeriv_for(Xtem,Ytem,K,sgx,sgy,EPS,Wt);%,NDIV);
        end
    BMtotal=sum(BM,3);
    [B_out,LL]=eigs(BMtotal,K);
%     [B2]=KernelDeriv(X,Y,K,sgx,sgy,EPS,Wt);
%     fprintf('sgx = %.3f  sgy = %.3f eps=%e\n', sgx,sgy,EPS);
    else    
        [B_out]=KernelDeriv_single(X,Y,K,sgx,sgy,EPS,SSdata,obsY);%,NDIV);

    end
end

% fprintf(Bm,'\nB_out matrix\n%3f\n',B);
% fclose(Bm);
% disp(B_out);
endtime = cputime - starttime;
%fprintf('Time = %3.4f\n', endtime);

end

