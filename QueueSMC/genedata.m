totaltime=cputime;

%%%
%random
a=clock;
seed=floor(a(6));
s = RandStream('mt19937ar','Seed',seed);
RandStream.setGlobalStream(s);
j=1;
while(j<=100)   
theta1temp = 10*rand(1,1);
if(theta1temp>2 && theta1temp<8)
    theta1(j)=theta1temp;
    j=j+1;
end
end

theta2=theta1+10*rand(1,100);
j=1;
while(j<=100)   
theta3t=1/3.*rand(1,1);
if(theta3t>0.1 && theta3t<0.3)
    theta3(j)=theta3t;
    j=j+1;
end
end

obsSS=zeros(100,22);
for i=1:100
    obsdata(i,:)=simulationMG1([theta1(i) theta2(i) theta3(i)]);
%     SumStatQuntSim = quantile(obsdata(i,:),20);
%     SumStatMaxSim = max(obsdata(i,:));
%     SumStatMinSim = min(obsdata(i,:));
%     obsSS(i,:) = [SumStatMaxSim SumStatMinSim SumStatQuntSim];
end

save obsdata1.mat obsdata theta1 theta2 theta3;
