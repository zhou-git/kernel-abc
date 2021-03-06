function [ out ] = para( input )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Prior1 = 1;                     %Parameter 'Prior1' #for ground data
Prior2 = 5;                     %Prior2
Prior3 =0.2;                    %prior3

%%%%%%%%%%%%%%%%%%%%%
num_para = 3;
 bandwidth =[20;15;10;8;6;5;3;1;0.5;0.3];
%bandwidth =[2 2 2;1.5 1.5 1.5;1.2 1.2 1.2;1 1 1;0.7 0.7 0.7;0.5 0.5 0.5;0.3 0.3 0.3];%s 12;%7.2;                  %Bandwidth for rejction methods
%bandwidth =[2;1.5;1;0.8;0.6;0.45;0.3];%theta1
% bandwidth =[2;1.7;1.3;1.1;1;0.9;0.7];%theta2
%bandwidth =[1;0.5;0.3;0.1;0.05;0.03;0.01];theta3

%%%%%%%%%%%%%%%%%%%%%
num_quantile = 10;
num_summarys = num_quantile + 2;

%%%%Simulation Modes
MCMC = 0;
REJ = 1;
LLR = 0;                        %Local linear Regression on/off
GKDR= 0;
SEMI = 0;

%%for GKDR
eff_dim=3;
training=20000;
local_band=[15];%;30;30];

switch input
    case 'pri1'
        out = Prior1;        
    case 'pri2'
        out = Prior2;        
    case 'pri3'
        out = Prior3;      
    case 'bandwidth'
        out = bandwidth.*1.5;%for rej methods        
    case 'sim_num'
        out = simulation_num;    
    case 'c'
        out = 0.75;%calcc(bandwidth);
    case 'interations'
        out = iteration;
    case 'np'
        out = num_para;
    case 'LLR'
        out = LLR;
    case 'n_summarystatics'
        out = num_summarys;
    case 'MC_iters'
        out = 2000;     %%for MCMC method, the iterations of the chain
    case 'sigma'
        out = [0.5;1;0.1];        %step length for MCMC methods.S
    case 'n_data_semiauto_ssEstimation'
        out = 10000;
    case 'num_quantile'
        out = num_quantile;
    case 'gkdr'
        out = GKDR;
    case 'REJ'
        out = REJ;
    case 'MCMC'
        out = MCMC;
    case 'SEMI'
        out = SEMI;
    case 'effective_dimension'
        out = eff_dim;
    case 'training_number'
        out = training;
    case 'local_b'
        out = local_band;
end

end

