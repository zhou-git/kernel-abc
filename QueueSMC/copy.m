% recopy trajectories according to resampling results
function [theta, x]= copy(thetas,xs,N_sons)

theta=zeros(size(thetas));
x=zeros(size(xs));
N=length(N_sons);
ind=1;

for i=1:N
   
   % if copy then keep it here
   if (N_sons(1,i)>0)
      
      aux_x=xs(i,:,:);
      aux_theta=thetas(:,i);
 
      for j=ind:ind+N_sons(1,i)-1
         
         x(j,:,:)=aux_x;
         theta(:,j)=aux_theta;
         
      end
      
   end
   
   ind=ind+N_sons(1,i);
   
end





