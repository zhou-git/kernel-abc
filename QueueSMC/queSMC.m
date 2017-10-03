
%Last Modified: 2017/4/26

%ABC  Queue model main program

%F='rej_iter50';
%file = fopen(F,'a'); 
%profile on

tic1=cputime;
plot_switch=0;
%%parameter
P=12; 
M=5;               % number of pseudo-obs for each particle
alpha=0.95;         % percentage particles surviving
np=3;
N=100000;             % number of particles
epsilontarget=0.3; % final target
    
C = {'k','b','r','g','y'};
D = {'o','+','*','.','x'};
% theta=zeros(1,N);   % particles
a=clock;
seed=floor(a(6));
s = RandStream('mt19937ar','Seed',seed);
RandStream.setGlobalStream(s);


fprintf('\n\nQueue Model...SMC\n\n epsilon %f, N %d',epsilontarget,N);
load('obsdata.mat');            %theta1, th2,th3 and obsdata obsSS
obs1 = quantile(obsdata,P-2,2);
obs2 = max(obsdata,[],2);                 %sampdata:Nx50
obs3 = min(obsdata,[],2);
obsSS=[obs2 obs3 obs1];
%load('Y.mat');
q=zeros(10,3);                      %acceptance rate
err=zeros(10,3);                    %msq error
mean=zeros(10,3);                   %empirical mean
%n_reps=300000;                      %number of simulation ,should change in smc
pribd=[10 20 0.333];                %boundness of prior and proposals

%Rssei=zeros(3,Acc);
round=10;
P=12;                               %number of summary statistics
x_dis=zeros(M,N);
Rssei=zeros(np,round);
prop_step=[0.5 0.5 1];
%ABC main process


for i=1:round
    obstemp=obsSS(i,:);
    obstheta=[theta1(i); theta2(i); theta3(i)]; %%prepare observed data: 1x3
    fprintf('\n True parameter: %f %f %f\n',obstheta);
    %%here begin the smc
    tic2=cputime;
   %%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%% Buffer setup for reproducibility
    counter=0;          % count the number of resampling
    epsilonstore=0;     % store epsilons
    essstore=0;
    accstore=0;
    esstresh=N/2;
    num_par_store=0;
    theta_resample=zeros(np,N,10);
    %%%%%%%%%%%%%%%%%%%%%%%%%%

    acc=zeros(1,N);     % acceptance rate MH
   %w=zeros(1,N);       % importance weights
    npa=zeros(1,N);     % active or not pseudo-obs for each particle
    pa=zeros(1,N);      % particle active or not, active if npa>0
    x=zeros(N,P,M);       % pseudo-observations

    theta_old=zeros(3,N); % old particles
    theta_prop=zeros(3,N); % candidate particles
   % npa_old=zeros(1,N);
    npa_prop=zeros(1,N);
    x_old=zeros(N,P,M);
    x_prop=zeros(N,P,M);
    tpa_old=0;
    reflevel=0;

    d=zeros(N,P,M);       % distance to true observations

    epsilon=100;        % just has to be superior to final schedule
    epsilonold=0;
    std_prop=0;         % std of uniform proposals

    %%initialization   sample particles
    pri1 = 10*rand(1,N);                        %theta:3xN
    pri2=pri1+10*rand(1,N);
    pri3=1/3.*rand(1,N);
    theta=[pri1; pri2; pri3];                   %3xN
    for i1=1:M                                  %M different f(x|theta),theta fixed
        sampdata=simulationMG1(theta,N);                    %main sampler 
        SumStatQuntSim = quantile(sampdata,P-2,2);
        SumStatMaxSim = max(sampdata,[],2);                 %sampdata:Nx50
        SumStatMinSim = min(sampdata,[],2);
        x(:,:,i1) = [SumStatMaxSim SumStatMinSim SumStatQuntSim];
    end
    x_old=x;
    SS_m=repmat(obstemp,[N 1 M]);
    for i2=1:M
        x_dis(i2,:)=rssq(x(:,:,i2)-SS_m(:,:,i2),2);  %distance between x and obsSS, root of squared sum
    end
    w=ones(1,N)./N;
    npa_old=M.*ones(1,N); 

    counter2=0;
    counter4=0;
    c3=0;
    
    stack_flag=0;
    stack_c=0;
    
    while (epsilon>epsilontarget)
        counter2=counter2+1;
        tic3=cputime;
        % find next level
        epsilonold=epsilon;
        reflevel=alpha*tpa(epsilonold,x_dis);
        
        fun = @(epsilon) tpa(epsilon,x_dis)-reflevel;
        epsilon=fzero(fun,[0 epsilonold]);
        if (epsilon<epsilontarget)
            epsilon=epsilontarget;
        end
        epsilonstore=[epsilonstore epsilon];
        essstore=[essstore 1/sum(w.^2)];
        [~,num_pa_tmp]=size(unique(theta(1,:)));
        num_par_store=[num_par_store num_pa_tmp];

        %%compute the associated weights and resample if necessary
        d=abs(x_dis);
        npa_old=sum((d<epsilonold),1);
        npa=sum((d<epsilon),1);
        a=find(npa_old>0); 
        b=find(npa_old==0);
        w(1,a)=w(1,a).*npa(1,a)./npa_old(1,a);
        w(1,b)=zeros(1,length(b));
        w=w./sum(w);

        %%resample if necessary and stroe the particals 
        %if (tpa(epsilon,x)<beta)
        if (sum(w.^2)*esstresh>1)
            N_sons=rsdet(w);
            [theta,x]=copy(theta,x,N_sons);  %%size(unique(theta(1,:)))
            w=ones(1,N)./N;
            counter=counter+1;
            Rssei_tmp = minus(sum(theta,2)./N,obstheta).^2;
            fprintf('\n epsilon %f', epsilon);
            fprintf(' rssei %f %f %f',Rssei_tmp);
            if (mod(counter,5)==0 && plot_switch==1 && i==1)
                for k = 1:3 
                    figure(k); 
                    histfit(theta(k,:),20,'kernel'); 
                    temp=['semitheta ',num2str(k),' at resample: ',num2str(counter),'.eps'];                    
                    temp2=['semitheta ',num2str(k),' at resample: ',num2str(counter)];      
                    xlabel(temp2);
                    saveas(gcf,temp);
                end 
            end
        end

        %%move particles around using MCMC steps
        theta_old=theta;                            % recopy values 
        x_dis_old=x_dis;
        x_old=x;

        %%propose for each particle but could only propose for alive particle to avoid wasting time
        d_old=abs(x_dis_old);
        npa_old=sum((d_old<epsilon),1); % compute number of active pseudo-obs for each particle
        indexalive=find(npa_old>0);

        %%adjust std of proposal and sample proposal for MH step
        %std_prop = zeros(3);
        th_prop=[0 10;10 20;0 0.4];
        for i3=[1 3]
            std_prop=prop_step(i3)*sqrt(sum(w.*(theta(i3,:).^2))-sum(w.*theta(i3,:))^2);
            theta_prop(i3,:)=theta_old(i3,:)+std_prop.*randn(1,N);  % sample candidate for MH steps
        end
        theta_prop = abs(theta_prop);
        std_prop=prop_step(2)*sqrt(sum(w.*((theta(2,:)-theta(1,:)).^2))-sum(w.*(theta(2,:)-theta(1,:)))^2);
        theta2_prop= theta(2,:)-theta(1,:)+std_prop.*randn(1,N);
        theta_prop(2,:)=theta_prop(1,:)+theta2_prop;
        %std_prop=0.15;
        
        
        for i4=1:M                                  %M different f(x|theta),theta fixed
            sampdata=simulationMG1(theta_prop,N);                    %main sampler 
            SumStatQuntSim = quantile(sampdata,P-2,2);
            SumStatMaxSim = max(sampdata,[],2);                 %sampdata:Nx50
            SumStatMinSim = min(sampdata,[],2);
            x_prop(:,:,i4) = [SumStatMaxSim SumStatMinSim SumStatQuntSim];
        end
        SS_m=repmat(obstemp,[N 1 M]);
        for i2=1:M
            x_dis(i2,:)=rssq(x_prop(:,:,i2)-SS_m(:,:,i2),2);  %distance between x and obsSS, root of squared sum
        end
        %w=ones(1,N)./N;
        %npa_old=M.*ones(1,N); 
        %d=abs(x_prop);
        npa_prop=sum((x_dis<epsilon),1);
        if(npa_prop==0)
            stack_c=stack_c+1;
            if(stack_c==2000)
                stack_flag=0;
                fprintf('break!\n');
                break;
            end
        end

        %%compute MH acceptance ratio
        acc=zeros(1,N);
        acc(1,indexalive)=npa_prop(1,indexalive)./npa_old(1,indexalive);   % acceptance proba for MH
        aq=find(abs(theta_prop(1,:)>pribd(1) | abs(theta_prop(2,:))>pribd(2)|abs(theta_prop(3,:))>pribd(3)));
        aq=unique(aq);
        acc(1,aq)=0;%zeros(1,length(aq));
        u=rand(1,N);
        aq=find(u<acc);
        accstore=[accstore length(aq)/length(indexalive)];
        accrate_store=length(aq)/length(indexalive);

        %%update accepted proposals
        theta(:,aq)=theta_prop(:,aq);                 % modify the values for accepted MH prop
        x(aq,:,:)=x_prop(aq,:,:);
        
        for i2=1:M
            x_dis(i2,:)=rssq(x(:,:,i2)-SS_m(:,:,i2),2);  %distance between x and obsSS, root of squared sum
        end
        d=abs(x_dis);
        
        if(accrate_store <= 0.0002)
%            break;
        end
    end
    % resample particles, build normalized histogram and target pdf
    N_sons=rsdet(w);
    [theta,x]=copy(theta,x,N_sons);
    if(i==1&&plot_switch==1)
        for k = 1:3 
            figure(k); 
            histfit(theta(k,:),20,'kernel'); 
            temp=['rejtheta ',num2str(k),'.eps'];            
            temp2=['rejtheta ',num2str(k)];      
            xlabel(temp2);
            saveas(gcf,temp);
        end
        % display epsilon and ess
        figure(4)
        plot(epsilonstore(1,2:end));
        saveas(gcf,'rejepsilon.eps');
        figure(5)
        plot(essstore(1,2:end));
        saveas(gcf,'rejESS.eps');
        figure(6)
        plot(accstore(1,2:end));
        saveas(gcf,'rejAcceptance rate.eps');
        figure(7)
        plot(num_par_store(1,2:end));
        saveas(gcf,'rejAcceptance rate.eps');
    end
    toc2=cputime-tic2;
    fprintf('\nRound: %d iteration: %d cpu time:  %d epsilon: %f Unique partical: %d\n',i, counter2,toc2,epsilon,num_pa_tmp);
    Rssei(:,i) = minus(sum(theta,2)./N,obstheta).^2;
    fprintf('Rssei %f %f %f\n****************************\n',Rssei(:,i));
end
RSSE=sum(Rssei,2)./round;
rssestd=std(Rssei,0,2);
tot_time=cputime-tic1;
fprintf('\n Total time spent is %d\n',tot_time);
fprintf('#############################\n\nFinal RSSE: %f %f %f\n STD: %f %f %f\n',RSSE,rssestd);

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
