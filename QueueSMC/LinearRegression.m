
function [betas] = LinearRegression(theta_buffer,s_buffer,i,ns,power,np,epi)
    t_f = zeros(i,np);
    s_f = zeros(ns*power,i);
    N=ns+1;
    for j=1:i
        t_f(j,:) = theta_buffer(j,:); 
        s_f(:,j) = s_buffer(:,j);
    end
    x = horzcat(ones(i,1),s_f');
    betas = (x'*x+epi*N.*eye(N))\(x'*t_f);
    Beta = betas(2:N,:);
end


