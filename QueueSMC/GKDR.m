function [B_out,V,para_gkdr] = GKDR(G,fix,obstemp,focuspara,band_fix)
%GKDR Summary of this function goes here
%   Detailed explanation goes here
% Bm=fopen('Bmatrix_s.txt','a+');
Np =20000;
Nt=20000;
N=1000; %Np>=4*N according to cv
np=3;
weight_k=1;
c=0.75;
P=size(obstemp,2);
CV_num=1000;
%% generate training data
pri1 = 10*rand(1,Np);                        %theta:3xN
pri2=pri1+10*rand(1,Np);
pri3=1/3.*rand(1,Np);
theta=[pri1' pri2' pri3'];
sampdata=simulationMG1(theta',Np);                    %main sampler 
SumStatQuntSim = quantile(sampdata,P-2,2);
SumStatMaxSim = max(sampdata,[],2);                 %sampdata:Nx50
SumStatMinSim = min(sampdata,[],2);
x = [SumStatMaxSim SumStatMinSim SumStatQuntSim];
tsyE=x';
%tsyE=abs((tsyE-mean(tsyE))./std(tsyE));
% test data
pri1t = 10*rand(1,Nt);                        %theta:3xN
pri2t=pri1t+10*rand(1,Nt);
pri3t=1/3.*rand(1,Nt);
thetat=[pri1t' pri2t' pri3t'];
sampdatat=simulationMG1(thetat',Nt);                    %main sampler 
SumStatQuntSimt = quantile(sampdatat,P-2,2);
SumStatMaxSimt = max(sampdatat,[],2);                 %sampdata:Nx50
SumStatMinSimt = min(sampdatat,[],2);
xt = [SumStatMaxSimt SumStatMinSimt SumStatQuntSimt];
tsyEt=xt';
%tsyEt=abs((tsyEt-mean(tsyEt))./std(tsyEt));
SumStatDat=obstemp;
K = G;

thetas=theta;
[M,tnum]=size(tsyE);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%Simulation     1,no rejection 2.with rejection
%N: number of training data for all parameters(different B)
%1:
starttime=cputime;

%%%calculate distance and random permute
% 
% rtnum=randperm(tnum);
% thetas=thetas(rtnum,:);
% tsyE=tsyE(:,rtnum);

if(fix==0)
    fprintf('Cross validation ...\n');
    fprintf('#sample = %d, dim of X = %d, effective dim = %d\n', N, M,G);
    candx=[0.25 0.5 0.75 1 2];  % candidates for CV
    candy = [0.25 0.5 1 2];
    eps = [0.0001 0.001];%[0.00001];

    % For cross-validation
%     NCV=5;      % Number of cross-validation
%     DEC=10;
%     lx=ceil(N/NCV);
%     ei=cumsum(lx.*ones(1,NCV),2);  
%     si=ei-(lx-1).*ones(1,NCV);
%     ei(NCV)=N;       % si: staring idx, ei: ending idx
    err_tbl=zeros(length(candx)*length(candy)*length(eps), 1);
    rsse_b=zeros(1,Np);

    for in=1:Np                 %(round(N/size_data))
        SumStatSim = tsyE(:,in);    %[SumStatMaxSim SumStatMinSim SumStatQuntSim]';
        rsse_b(in) = minus(SumStatSim,SumStatDat')'*minus(SumStatSim,SumStatDat'); 
    end

    [rsse_t,I]=sort(rsse_b);
    Ym=thetas(I,:);
    %for cv
   % Y_avr=sum(Ym((1:1000),:),1)./1000;
    Xm=tsyE(:,I)';
    Ytt=Ym(1:N,:);
%     Xtt=Xm(1:N,:);
    X=Xm(1:N,:);
    rtt=rsse_t(1:N);
    rsseup=rtt(N);
    W1=1-(rtt/rsseup);
    W1(W1==0)=0.00001;
    W2 = c.*(1.-(rtt./rsseup).^2);
    W3=(1-(rtt./rsseup).^2).^3;
%     Xstd=std(Xtt);
%     Xstd(Xstd==0)=1;
%    Ystd=std(Ytt);
%     X=Xtt./repmat(Xstd,N,1);
%     %Yttm=mean(Y2);
    Y_t=Ytt;%abs((Ytt-mean(Ytt))./std(Ytt));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if(weight_k==1)
        W=W1';
    elseif(weight_k==2)
        W=W2';
    elseif(weight_k==3)
        W=W3';
    elseif(weight_k==0)
        W=ones(1,Acc);
    end
    if(focuspara==0)
        Y=Y_t;
    elseif(focuspara==1)
        Y=Y_t(:,1);
    elseif(focuspara==2)
        Y=Y_t(:,2);
    elseif(focuspara==3)
        Y=Y_t(:,3);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    X=(X-repmat(mean(X),N,1));
    Xstd=std(X);
    Xstd(Xstd==0)=1;
    X=abs(X./repmat(Xstd,N,1));
    Y=(Y-repmat(mean(Y),N,1));
    Ystd=std(Y);
    Y=abs(Y./repmat(Ystd,N,1));

    sgx0=MedianDist(X);   % Basic value for bandwidth for X
    sgy0=MedianDist(Y);   % Basic value for bandwidth for Y 
    Ystd=std(Ytt);
    ridx=randperm(N);  % random order 
    Xr=X(ridx,:);
    Yr=Y(ridx,:);  
    Wr=W(ridx,:); 
    %% parallel validation
    sgx_v=sgx0.*candx;
    sgy_v=sgy0.*candy;
    tt1=length(candx)*length(candy)*length(eps);
    tt2=length(candy)*length(eps);
    for hhh=1:tt1
        opth1=ceil(hhh/tt2);
        sgx=sgx_v(opth1);
        rr1=hhh-(opth1-1)*tt2;
        optk1=ceil(rr1/length(eps));
        sgy=sgy_v(optk1);
        opte1=mod(rr1,length(eps));
        if opte1==0
            opte1=1;
        end 
%     for h=1:length(candx)
%         sgx=sgx0*candx(h);% five different values for CV
%         for k=1:length(candy)
%           sgy=sgy0*candy(k);            
%           for ll=1:length(eps)
            EPS = eps(opte1);
%                 for i=1:NCV
%                     ri=si(i):ei(i);
%                     Xe=Xr; Ye=Yr; We=Wr;
%                     Xe(ri,:)=[];
%                     Ye(ri,:)=[];    % Xe, Ye: trainig sample for CV
%                     We(ri,:)=[];
%                     Xt=Xr(ri,:);
%                     Yt=Yr(ri,:);    % Xt, Yt: test sample for CV
            [B,t,V]=KernelDeriv_for(Xr,Yr,K,sgx,sgy,EPS,Wr);
            
% 
%                     sampdata=simulationMG1(Y_avr',1);                    %main sampler 
%                     SumStatQuntSim_tem = quantile(sampdata,P-2,2);
%                     SumStatMaxSim_tem = max(sampdata,[],2);                 %sampdata:Nx50
%                     SumStatMinSim_tem = min(sampdata,[],2);
%                     x_tem = [SumStatMaxSim_tem SumStatMinSim_tem SumStatQuntSim_tem];

              %  for in=1:10000                 %(round(N/size_data))
%                     SumStatSim = tsyE(1:10000,:)*B;    %[SumStatMaxSim SumStatMinSim SumStatQuntSim]';
%                     x_b=repmat(x_tem,10000,1)*B;
%                     rsse_tmp = minus(SumStatSim,x_b); 
%                     rsse_b1=rsst_tmp'*A*rsse_tmp;
%                 %    end
%                     [~,I]=sort(rsse_b1);
%                     Ym1=thetas(I,:);
%                     Ym2=sum(Ym1(1:1000,:),1)./1000;
%                     dd=abs(Y_avr-Ym2);
            % kNN regression for CV
            Yo=zeros(N,np);
            nnidx=knnsearch(tsyEt'*B,Xr*B, 'K', CV_num, 'NSMethod', 'kdtree'); 
            for jjj=1:N
                indx=nnidx(jjj,:);
                Yo(jjj,:)=sum(thetat(indx',:),1)./CV_num;%Ye
            end
            dd=Yr-Yo;
         %  inx_err=(h-1)*length(candy)*length(eps)+(k-1)*length(eps)+ll;
            err_tbl(hhh)=sum(sum(dd.*dd,1),2);                    
            %5 means divide the data to 5 parts. 
%                 end

%           end
%         end
    end

    %     fprintf('Validation Done\n');

    [c, midx]=min(err_tbl);
    opth=ceil(midx/(length(candy)*length(eps)));
    rr=midx-(opth-1)*length(candy)*length(eps);
    optk=ceil(rr/length(eps));
    opte=mod(rr,length(eps));
    if opte==0
        opte=length(eps);
    end   
    
    %save errortable.mat err_tbl;
%%%%%%%%%%%%%%%%%%%%%%%%
%set up X,Y and weight again as opimum value
    
    rsse_b=zeros(1,Np);
    for in=1:Np  %(round(N/size_data))
        SumStatSim = tsyE(:,in);%[SumStatMaxSim SumStatMinSim SumStatQuntSim]';
        rsse_b(in) = minus(SumStatSim,SumStatDat')'*minus(SumStatSim,SumStatDat'); 
    end
    [rsse_t,I]=sort(rsse_b);
    Ym=thetas(I,:);
    Xm=tsyE(:,I)';
    Ytt=Ym(1:N,:);
    Xtt=Xm(1:N,:);
    rtt=rsse_t(1:N);
    rsseup=rtt(N);
    W1=1-(rtt/rsseup);
    W1(W1==0)=0.00001;
    W2 = c.*(1.-(rtt./rsseup).^2);
    W3=(1-(rtt./rsseup).^2).^3;
    
    Xtt=(Xtt-repmat(mean(Xtt),N,1));
    Xstd=std(Xtt);
    Xstd(Xstd==0)=1;
    X=abs(Xtt./repmat(Xstd,N,1));
    
    Ytt=(Ytt-repmat(mean(Ytt),N,1));
    Ystd=std(Ytt);
    Ystd(Ystd==0)=1;
    Ytt=abs(Ytt./repmat(Ystd,N,1));
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if(weight_k==1)
        W=W1';
    elseif(weight_k==2)
        W=W2';
    elseif(weight_k==3)
        W=W3';
    elseif(weight_k==0)
        W=ones(1,Acc);
    end
    if(focuspara==0)
        Y=Ytt;
    elseif(focuspara==1)
        Y=Ytt(:,1);
    elseif(focuspara==2)
        Y=Ytt(:,2);
    elseif(focuspara==3)
        Y=Ytt(:,3);
    end    
    %%%%%%%%%%%%%%%%%%%%%%%
    sgx=sgx0*candx(opth);
    sgy=sgy0*candy(optk);
    EPS=eps(opte);
    [B,~,V]=KernelDeriv(X,Y,K,sgx,sgy,EPS,W);
    para_gkdr=[sgx sgy EPS];
    B_out=B;
    %   fprintf('NCV=%d, DEC=%d, K=%d\n', NCV, DEC);
    fprintf(' Npop: %d,N: %d',Np,N);
    fprintf(',G= %d sgx=%f,sgy=%f,eps:%f,W:%d,Yfocus:%d\n',K,sgx,sgy,EPS,weight_k,focuspara);
    
%     fprintf('Extimation Error: %f\n', mean(err));
else
%     Y = Y_t(:,1);
    Npcand=Np;
    rsse_b=zeros(1,Npcand);
    for in=1:Npcand  %(round(N/size_data))
        SumStatSim = tsyE(:,in);%[SumStatMaxSim SumStatMinSim SumStatQuntSim]';
        rsse_b(in) = minus(SumStatSim,SumStatDat')'*minus(SumStatSim,SumStatDat'); 
    end
    [rsse_t,I]=sort(rsse_b);
    Ym=thetas(I,:);
    Xm=tsyE(:,I)';
    Ytt=Ym(1:N,:);
    Xtt=Xm(1:N,:);
    rtt=rsse_t(1:N);
    rsseup=rtt(N);
    W1=1-(rtt/rsseup);
    W1(W1==0)=0.00001;
    W2 = c.*(1.-(rtt./rsseup).^2);
    W3=(1-(rtt./rsseup).^2).^3;
    Xtt=(Xtt-repmat(mean(Xtt),N,1));
    Xstd=std(Xtt);
    Xstd(Xstd==0)=1;
    X=abs(Xtt./repmat(Xstd,N,1));
    
    Ytt=(Ytt-repmat(mean(Ytt),N,1));
    Ystd=std(Ytt);
    Ystd(Ystd==0)=1;
    Ytt=abs(Ytt./repmat(Ystd,N,1));
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if(weight_k==1)
        W=W1';
    elseif(weight_k==2)
        W=W2';
    elseif(weight_k==3)
        W=W3';
    elseif(weight_k==0)
        W=ones(1,Acc);
    end
    if(focuspara==0)
        Y=Ytt;
    elseif(focuspara==1)
        Y=Ytt(:,1);
    elseif(focuspara==2)
        Y=Ytt(:,2);
    elseif(focuspara==3)
        Y=Ytt(:,3);
    end    
    %%%%%%%%%%%%%%%%%%%%%%%
    sgx=band_fix(1);
    sgy=band_fix(2);
    EPS=band_fix(3);
    [B,~,V]=KernelDeriv(X,Y,K,sgx,sgy,EPS,W);
    
    B_out=B;
    para_gkdr=band_fix;
    
    %   fprintf('NCV=%d, DEC=%d, K=%d\n', NCV, DEC);
    fprintf(' Npop: %d,N: %d',Npcand,N);
    fprintf(',G= %d sgx=%f,sgy=%f,eps:%f,W:%d,Yfocus:%d\n',K,sgx,sgy,EPS,weight_k,focuspara);
end

% fprintf(Bm,'\nB_out matrix\n%3f\n',B);
% fclose(Bm);
% disp(B_out);
endtime = cputime - starttime;
fprintf('Time = %3.4f\n', endtime);

end

