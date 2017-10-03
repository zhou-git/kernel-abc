
function [thetam] = LocalLinearRgs(theta_buffer,s_buffer,w,i,ns)
    t_f = zeros(i,3);
    s_f = zeros(ns,i);
%    w_array_f = zeros(i,1);
    w_f = zeros(i);
    for j=1:i
        t_f(j,:) = theta_buffer(j,:); 
        s_f(:,j) = s_buffer(:,j);
%        w_array_f(j) = w(j);
        w_f(j,j) = w(j);
    end
    x = horzcat(ones(i,1),s_f');
    betas = inv(x'*w_f*x)*x'*w_f*t_f;
    thetam = betas(1,:)';
%     x = horzcat(ones(i,1),s_buffer');
%     betas = inv(x'*w*x)*x'*w*theta_buffer;
%     thetam = betas(1,:)';
end


