function [q_out,mser,mean,i] = RejGKDR(simnum,theta,c,delta,msdata,Y,obsdata,data_partation,size_data,B)
%Rejection method
%
Acc = 5000;

LLR = para('LLR');
G=para('effective_dimension');
N=round(simnum/size_data);
%%Buffer
diff_SS_buffer = zeros(G,Acc);    %%Memory run out, do not need so large matrix for weights
weight_buffer = zeros(Acc,1);
theta_buffer = zeros(Acc,1);
SS_buffer = zeros(G,Acc);
SS_bufferall = zeros(G,N);

SumStatDat = obsdata;
A=diag([(5.4/9.45)^(-1) (7.5/9.45)^(-1) (2/9.45)^(-1)]);
A=diag([(47/69)^(-1) (51/69)^(-1)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%
%Simulation

%load('ms1data_0_50_1.mat');%stheta
i=0;
j=0;
ptr=round((data_partation-1)*N/size_data);
while i<Acc && j<N
    j=j+1;
    SumStatSim = msdata(ptr+j,:);%[SumStatMaxSim SumStatMinSim SumStatQuntSim]';
    n1 =minus(SumStatSim*B,SumStatDat*B)*A*minus(SumStatSim*B,SumStatDat*B)'; %norm(minus(SumStatSim*B,SumStatDat*B));
    if(n1 <= delta)
        i = i + 1;
        theta_buffer(i) = Y(j+ptr);
        diff_SS_buffer(:,i) = SumStatSim*B-SumStatDat*B;
        SS_buffer(:,i)=SumStatSim*B;
        weight_buffer(i,1) = c/delta*(1-(n1/delta)^2);
    end
    SS_bufferall(:,j)=SumStatSim*B;
end
SS_bufferall=SS_bufferall(:,any(SS_bufferall));
if(i<Acc)    
    fprintf('number of accepted samples=%3d\n',i);
end
if (data_partation==1)
    fprintf('SS %f\n',sum(abs(SS_buffer),2)./i);
    fprintf('SSacc var %f\n', var(abs(SS_buffer),0,2));
    fprintf('SSall var %f\n', var(abs(SS_bufferall),0,2));
end
if(LLR)
    thetam = LocalLinearRgs(theta_buffer, diff_SS_buffer,weight_buffer,i,G);
    mser = ((thetam - theta).^2)./(theta.^2); 
    q = i/j;
else
    thetam = sum(theta_buffer);
    mean= thetam'/i;              %%use mean as point estimation
    mser = ((mean - theta).^2)./(theta.^2);
%     fprintf('stderr=%f\n',std(theta_buffer,0,1));
    q_out = i/j;
end 

end
