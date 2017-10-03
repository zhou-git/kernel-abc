function [ Beta ] = pilotrej(SSdata)
%PILOTMCMC Summary of this function goes here
%summary statistics are sorted interdeparture time
c=0.75;
% load('obsdata.mat');
% load('train.mat');
% 
% % X = zeros(N,M);   %t:N,ns:M
% % Y_t=zeros(N,np);
% 
% theta=tPri;
% np=size(theta,2);
power=1;
epi=0.00001;
N=30000;
Acc=3000;
[~,P]=size(SSdata);
np=3;
%% generate training data
pri1 = 10*rand(1,N);                        %theta:3xN
pri2=pri1+10*rand(1,N);
pri3=1/3.*rand(1,N);
theta=[pri1' pri2' pri3'];
sampdata=simulationMG1(theta',N);                    %main sampler 
SumStatQuntSim = quantile(sampdata,P-2,2);
SumStatMaxSim = max(sampdata,[],2);                 %sampdata:Nx50
SumStatMinSim = min(sampdata,[],2);
x = [SumStatMaxSim SumStatMinSim SumStatQuntSim];
tsyE=x';

[nd,nreps]=size(tsyE);
W0=ones(nreps,1);

%A=diag(1./var(tsyE,0,2));%% distance is squared
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
%Simulation     1,no rejection 2.with rejection
%N: number of training data for all parameters(different B)
%1:
starttime=cputime;
SumStatDat=SSdata;

% tnum=60000;
rsse_b=zeros(1,nreps);
rsse_t=zeros(1,nreps);


for in=1:nreps%(round(N/size_data))

    SumStatSim = tsyE(:,in);%[SumStatMaxSim SumStatMinSim SumStatQuntSim]';
    rsse_b(in) = minus(SumStatSim,SumStatDat')'*minus(SumStatSim,SumStatDat'); 
end
stdE=var(tsyE,0,2);

%stdEm=repmat(stdE,1,Acc);
[rssY,order] = sort(rsse_b);
% display(order);
%temp2=tsyE(:,order);
temp2=x(order,:);

%xball=normrnd(0,1,nd,Acc);
%normmat=repmat(sqrt(sum((xball./stdEm).^2)),nd,1);
%xsamt=xball./normmat;
%xsam=xsamt.*repmat(rssY(Acc),nd,Acc);
%xsam=xsamt2.*stdE2m.*repmat(temp2(:,Acc),1,Acc);
% sse=sum(xsam.^2);
% idx=find(sse>rssY(Acc));

temp3=sort(temp2(1:Acc,:),2);%+xsam;
%tsyt=(temp3-tsyE)./stdEm;
tsyt=(temp3');%./stdEm;

thetaYt = theta(order,:);
thetaYacc = thetaYt(1:Acc,:);
tsyt =tsyt-repmat(mean(tsyt,2),1,Acc);
std1=std(tsyt,0,2);
tsyt=abs(tsyt./repmat(std1,1,Acc));
thetaYacc = thetaYacc-repmat(mean(thetaYacc),Acc,1);
std2=std(thetaYacc);
thetaYacc =abs(thetaYacc./repmat(std2,Acc,1));

    %Calculate SS
    betas=LinearRegression(thetaYacc,tsyt,Acc,P,power,np,epi);
%   betas= LocalLinearRgs(para_buff, SS_buffer,weight_buffer,ss,ns,np);
    Beta = betas(2:P*power+1,:);
end

