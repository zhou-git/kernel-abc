function [ data ] = simulationMG1(state)
%SIMULATIONMG1 Summary of this function goes here
%simulation of M/G/1 quene with three parameters.theta1, theta2, theta3.
%service time is uniformly distributed in[theta1,theta2].
%interarrival time is expontional distributed with theta3.
%   Detailed explanation goes here
%data: row vector
    
%%

data = zeros(1,50);
theta1 =state(1);
theta2 =state(2);
theta3 = state(3);

U = theta1 + (theta2 - theta1)*rand(1,1);
YTotal = U;
data(1) = U;
IntArrTotal = -log(rand(1))/theta3;


for i=2:50;
    U = theta1 + (theta2 - theta1)*rand(1,1);
    W = -log(rand(1))/theta3;    
    IntArrTotal = IntArrTotal + W;

    if(IntArrTotal <= YTotal)  
        Y = U;
    else
        Y = U + IntArrTotal - YTotal;
    end
    
    YTotal = YTotal + Y;
     data(i) = Y;
end

end


