function [ Yerr ] = sortsampGKDR( theta,obsdata,msdata,Y,B )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Acc = para('Acc');%3000;
N=100000;
M=7;
np=para('np');
LLR = para('LLR');
G=para('n_summarystatics');

madv=mad(msdata*B);
A=diag(1./madv);
rsse_b=zeros(1,N);
rsse_t=zeros(1,N);
%%%%%%%%%%%%%%%%%%%%%%%%%
%Prepare Data

SumStatDat = obsdata*B;

%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:N;%(round(N/size_data))

    SumStatSim = msdata(i,:)*B;
    rsse_b(i) = minus(SumStatSim,SumStatDat)*minus(SumStatSim,SumStatDat)'; %

end
[rsse_t,I]=sort(rsse_b);
Yt=Y(I);
Yerr=norm(minus(mean(Yt(1:Acc)),theta));

end

