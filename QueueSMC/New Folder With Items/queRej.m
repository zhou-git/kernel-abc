%Last Modified: 31 August.
%ABC  Queue model main program

%F='rej_iter50';
%file = fopen(F,'a'); 
%profile on

totaltime=cputime;

%%%
%random
a=clock;
seed=floor(a(6));
s = RandStream('mt19937ar','Seed',seed);
RandStream.setGlobalStream(s);

Acc=3000;
fprintf('\n\nQueue Model...REJECTION\n\n');
load('obsdata.mat');
load('Y.mat');
q=zeros(10,3);
err=zeros(10,3);
mean=zeros(10,3);
n_reps=300000;
thetaY=[pri1' pri2' pri3'];
Rssei=zeros(3,Acc);
round=30;
%%%
%ABC main process
%%%
%1.Rejection method(GKDR, LLR)
for i=1:round
    obstemp=obsSS(i,:);
    obstheta=[theta1(i); theta2(i); theta3(i)];
    Ytemp=YSS;
    temp=minus(repmat(obstemp, n_reps,1), Ytemp).^2;
    se = abs(sum(temp,2));
    [rssY,order] = sort(se);
    thetaYt = thetaY(order,:);
    thetaYacc = thetaYt(1:Acc,:)';
    thetas=repmat(obstheta,1,Acc);
   % Rssei(:,i) =sum(minus(thetaYacc,thetas).^2,2);
    Rssei(:,i)=minus(sum(thetaYacc,2)./Acc,obstheta).^2;
    fprintf(' %f %f %f\n',Rssei(:,i));
end
RSSE=sum(Rssei,2)./round;
fprintf('\nRSSE is %f \n',RSSE);

%%%
totaltimeend=cputime-totaltime;
fprintf('Total time spent = %3f\n',totaltimeend);
fprintf('END of ABC\n\n');
% rmse = sum(rmse_tmp,2)./err_count;
% q_final = sum(qs)/ti;
% fprintf('%12s%3f %3f %3f\n', 'final rmse is:',rmse);
% fprintf('%10s %3f\n','q is:',q_final);
% fprintf('\nFinal rmse = %3f %3f %3f',rmse);
% fprintf('\nAcceptance rate is %3f',q_final);
% dlmwrite('rmse.txt',rmse);
% fprintf(file,'%3f  ',q_final);
%fclose(file);
% quit;
%profile viewer
%p = profile('info');
%profsave(p,'profile_results');
