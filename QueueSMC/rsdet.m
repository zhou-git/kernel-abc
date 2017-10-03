% deterministic selection procedure
% (Kitagawa JCGS 1996)
function N_sons= rsdet(weight)

N=length(weight);
u=zeros(1,N);
N_sons=zeros(1,N);

% generate the cumulative distribution
dist=cumsum(weight);
aux=rand(1);

u=aux:1:(N-1+aux);
u=u./N;
j=1;
for i=1:N
   while (u(1,i)>dist(1,j))
      j=j+1;
   end
   N_sons(1,j)=N_sons(1,j)+1;
end



