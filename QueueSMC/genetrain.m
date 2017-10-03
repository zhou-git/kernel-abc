totaltime=cputime;

ntotal=10000;
%%%
%random
a=clock;
seed=floor(a(6));
s = RandStream('mt19937ar','Seed',seed);
RandStream.setGlobalStream(s);

pri1 = 10*rand(1,ntotal);
pri2=pri1+10*rand(1,ntotal);
pri3=1/3.*rand(1,ntotal);
tPri=[pri1;pri2;pri3]';
tSS=zeros(ntotal,22);
tY=zeros(ntotal,50);
parfor i=1:ntotal
    tY(i,:)=simulationMG1([pri1(i) pri2(i) pri3(i)]);
    SumStatQuntSim = quantile(tY(i,:),20);
    SumStatMaxSim = max(tY(i,:));
    SumStatMinSim = min(tY(i,:));
    tSS(i,:) = [SumStatMaxSim SumStatMinSim SumStatQuntSim]';
end

save train.mat tY tPri tSS;
