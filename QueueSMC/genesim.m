totaltime=cputime;

ntotal=300000;
%%%
%random
a=clock;
seed=floor(a(6));
s = RandStream('mt19937ar','Seed',seed);
RandStream.setGlobalStream(s);

pri1 = 10*rand(1,ntotal);
pri2=pri1+10*rand(1,ntotal);
pri3=1/3.*rand(1,ntotal);
YSS=zeros(ntotal,22);
Y=zeros(ntotal,50);
parfor i=1:ntotal
    Y(i,:)=simulationMG1([pri1(i) pri2(i) pri3(i)]);
     SumStatQuntSim = quantile(Y(i,:),20);
    SumStatMaxSim = max(Y(i,:));
    SumStatMinSim = min(Y(i,:));
    YSS(i,:) = [SumStatMaxSim SumStatMinSim SumStatQuntSim]';
end

save Y.mat Y pri1 pri2 pri3 YSS;