%Last Modified: 31 August.
%ABC  MS model main program

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

%%%
%Parameters setup
n_t=50;
n_reps=1000000;
burnin=50;


type=1;%%  E0~E2
fprintf('\ntype %d',type);
fix=0;
%%%
%ABC main process
%%%
%1.Rejection method(GKDR, LLR)

load('sYE2_1.mat');
load('sYE2_2.mat');
load('sYE2_3.mat');

load('Y.mat');
load('yobsFix.mat');

% b1Y=cat(3,b1Y1,b1Y2);
sYE2=[sYE2_1;sYE2_2;sYE2_3];
thetaY=[log_r log_sig_e phi];
Rssei = zeros(3,100);
if(type==1)
    obs_syE=[obs_syE0];
    sYE=[sYE0];
elseif(type==2)
    obs_syE=obs_syE1;
    sYE=sYE1;
elseif(type==3)
    obs_syE=obs_y;
    sYE=Y;
elseif(type==4)
    obs_syE=[obs_syE1;obs_syE2(27:126,:);obs_syE2(375:424,:)];
    sYE=[sYE1;sYE2(27:126,:);sYE2(375:424,:)];
elseif(type==5)
    obs_syE=[obs_syE1;obs_y];
    sYE=[sYE1;Y];
end

% sYEm=mean(sYE,2);
% sYEabs=abs(sYE);
% sYEt=std(sYEabs,0,2);
% sYEf=(sYE-repmat(sYEm,1,n_reps))./repmat(sYEt,1,n_reps);
stdY=std(sYE,0,2);
stdY(stdY==0)=1;
% stdYmm=mean(sYE,2);
% stdYm=repmat(stdYmm,1,n_reps);
sYEf=sYE./repmat(stdY,1,n_reps);
for i = 1:30
    B = GKDR(5,fix,type,i);
    obs_syEf=obs_syE(:,i)./stdY;
    thetay = [obs_log_r(i) obs_log_sig_e(i) obs_phi(i)];
    obsyd=obs_syEf'*B;
    sYEd=sYEf'*B;
    temp=minus(repmat(obsyd, n_reps,1), sYEd).^2;
    se = abs(sum(temp,2));
    [rssY,order] = sort(se);
    thetaYt = thetaY(order,:);
    thetaYacc = thetaYt(1:100,:)';
    thetaMean = mean(thetaYacc,2);
    thetas=repmat(thetay',1,100);
   % epsilon=rssY(101);
%     krlwgt=3/4.*(1-(rssY(1:100)./epsilon).^2);
   % krlwgt=(35/32).*((1-(rssY(1:100)./epsilon).^2).^3);
   % kwmean=sum(krlwgt);
   % se=thetaYacc*krlwgt./kwmean;

   % Rssei(:,i) = abs(minus(se,thetay'));%mean(minus(thetaYacc,thetas).^2,2);
    Rssei(:,i) =sum(minus(thetaYacc,thetas).^2,2);

    fprintf('%f %f %f\n',Rssei(:,i));
end
RSSE=sum(Rssei,2)./30;
fprintf('RSSE is %f \n',RSSE);
%%%
totaltimeend=cputime-totaltime;
fprintf('Total time spent = %3f\n',totaltimeend);
fprintf('END of ABC\n\n');
