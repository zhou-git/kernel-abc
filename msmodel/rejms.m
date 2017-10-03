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
theta1 = para('pri1');
theta2 = para('pri2');
theta3 = para('pri3');
c = para('c');
np = para('np');
mcmc = para('MCMC');
rej = para('REJ');
semi = para('SEMI');
gkdr = para('gkdr');
llr=para('LLR');
theta=10;
fprintf('\n\nABC begin...\n\n');
%load('data.mat');
%%%
%ABC main process
%%%
%1.Rejection method(GKDR, LLR)
if(rej)
    delta_array = para('bandwidth');
    for hrej=1:5
        delta=delta_array(hrej);
        fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
        fprintf('Delta=%3.1f\n',delta);
        fprintf('Simulation Algorithm = REJECTION\n');
        fprintf('GKDR = %3.0f\n',gkdr);
        fprintf('LLR = %3.0f\n',llr);
        [q,err,mean] = sampms(theta,c,delta);  
        fprintf('MEAN = %3f %3f %3f\n',mean);
        fprintf('MSE = %3f %3f %3f\n',err);
        fprintf('Acceptance rate = %3f %3f %3f\n',q);
        fprintf('%%%%%%%%%%%END%%%%%%%%%%%%%%%\n\n');
    end
end

%%%
totaltimeend=cputime-totaltime;
fprintf('Total time spent = %3f\n',totaltimeend);
fprintf('END of ABC\n\n');