%Last Modified: 31 August.
%ABC  Queue model main program

%F='rej_iter50';
%file = fopen(F,'a'); 
%profile on
%clear all;
totaltime=cputime;

%%%
%random
a=clock;
seed=floor(a(6));
s = RandStream('mt19937ar','Seed',seed);
RandStream.setGlobalStream(s);

%%%
%Parameters setup
load('ms_20_1M_data.mat');
load('ms_20_1M_mad.mat');
load('ms_20_1M_Y.mat');
load('ms_20_1M_Ystd.mat');
load('Yobs.mat');
load('obsdata.mat');
np = para('np');
rej = para('REJ');
gkdr = para('gkdr');
llr=para('LLR');
fixed=1;
G=para('effective_dimension');
M=para('n_summarystatics');
tnum = para('training_number');
tdata_num=para('tdata_num');
index=1;
c=0.75;
ave=sum(ms_20_1M_data)./1000000;
fprintf('\n\nABC begin...\n\n');
fprintf('GKDR...\n');
fprintf('G: %f AccRate: %f \n',G,tnum/tdata_num);
mcmc=0;
err=zeros(1,100);
%ABC main process
%%%
%1.Rejection method(GKDR, LLR)
if(rej)
    if(fixed)
        load('ObAVE.mat');
        [B] = GKDR_fix(G,c,obsdata(1,:));
    else
        [B] = GKDR_abc(G,c,obsdata(1,:));
    end
    for pairs=1:100
        [err(pairs)]=sortsampGKDR(Yobs(pairs),obsdata(pairs,:),ms_20_1M_data,ms_20_1M_Y,B);
        fprintf('%f\n    ',err(pairs));
    end
%     for hrej=1:tt
%         delta=delta_array(hrej);
%         fprintf('Delta=%3.1f\n',delta);
%         for j=1:size_data
% %             [obsdata] = genedata(theta);
%             [q_b(j,:),~,mean_b(j,:),i_num(j)] = RejGKDR(simnum,theta,c,delta,msdata,Y,obsdataAve,j,size_data,B);  
%         end
%         mean_average=sum(mean_b)./size_data;
%         mean_stderr=std(mean_b);
%         q_average=sum(q_b)./size_data;
%         err_average=((mean_average - theta).^2)./(theta.^2);
%         Average_acceptance_num=sum(i_num)/size_data;
%          
%         fprintf('\nMEAN = %f, stderr=%f\n',mean_average,mean_stderr);
%         fprintf('MSE = %3.6f \n',err(pairs));
%         fprintf('AccRate = %3f, AccNum=%d\n\n',q_average(1),Average_acceptance_num);
%     end
    rsse=sum(err)/100;
    fprintf('\nRSSE is %f \n',rsse);
    fprintf('var %f',var(err));
end

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

% %%%
% if(semi)
%     for i = 1:ti
%         [q,err,stderr] = SemiAuto(data,theta);
%         
%     end
% end
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
