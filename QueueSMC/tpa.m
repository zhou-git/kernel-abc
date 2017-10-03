function t=tpa(epsilon,x)

N=size(x,2);
d=abs(x);
npa=sum((d<epsilon),1); % compute number of active pseudo-obs for each particle %%each row
pa=(npa>0);             % say whether particle active of not
t=sum(pa)/N;