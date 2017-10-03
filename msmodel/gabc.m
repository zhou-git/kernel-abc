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

%%%
%Parameters setup
theta1 = para('pri1');
theta2 = para('pri2');
theta3 = para('pri3');
np = para('np');
rej = para('REJ');
gkdr = para('gkdr');
llr=para('LLR');
fixed=1;
G=para('effective_dimension');
M=para('n_summarystatics');
theta = 10;
index=1;
c=0.75;
delta_array = para('bandwidth');
tt=size(delta_array,1);
err4=zeros(1,tt);
q4=zeros(1,tt);

fprintf('\n\nABC begin...\n\n');

%%%
%Simulated data
% data = simulationMG1(theta);

%%%
%ABC main process
%%%
%1.Rejection method(GKDR, LLR)
if(rej)
    if(fixed)
        [B] = GKDR_fix(G,c);
    else
        [B] = GKDR_abc(G,c);
    end
    for hrej=1:tt
        delta=delta_array(hrej);
        fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
        fprintf('Delta=%3.1f\n',delta);
        fprintf('Simulation Algorithm = REJECTION\n');
        fprintf('GKDR = %3.0f\n',gkdr);
        fprintf('LLR = %3.0f\n',llr);
        [q,err,mean] = RejGKDR(theta,c,delta,B);  
        fprintf('MEAN = %3f %3f %3f\n',mean);
        fprintf('MSE = %3f %3f %3f\n',err);
        fprintf('Acceptance rate = %3f %3f %3f\n',q);
        fprintf('%%%%%%%%%%%END%%%%%%%%%%%%%%%\n\n');
        err4(hrej)=err;
        q4(hrej)=q;
    end
end
%save err20_15 err1;
%save q20_15 q;
%%%
%2.MCMC method(GKDR,LLR)

if(mcmc)
    for ii=1:10
        fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
        fprintf('Data Size = %3.0f\n', sim_num*ii);
        fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
        fprintf('Simulation Algorithm = MCMC\n');
        fprintf('GKDR = %3.0f\n',gkdr);
        fprintf('LLR = %3.0f\n',llr);
        for i = 1:ti
            [q,err,mean] = McMc(data,theta,c);
            rmse_tmp(:,i) = err;
            mean_bf(:,i)=mean; 
            qs(i) = q;
        end
        stderr_final=std(mean_bf,0,2);
        rmse = sum(rmse_tmp,2)./ti;
        q_final = sum(qs)/ti;
        fprintf('MCMC OVERALL RESULTS...\n');
        fprintf('GKDR = %3.0f\n',gkdr);
        fprintf('LLR = %3.0f\n',llr);
        fprintf('Average MSE(50 iterations) = %3.4f %3.4f %3.4f\n',rmse);
        fprintf('STDERR = theta1=%3.4f,theta2=%3.4f,theta3=%3.4f\n', stderr_final);
        fprintf('Simulation Number for each iteration = %3.0f',sim_num*ii);
        fprintf('\nAcceptance rate is %3.4f\n',q_final);
    end
end

%%%
if(semi)
    for i = 1:ti
        [q,err,stderr] = SemiAuto(data,theta);
        
    end
end
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
