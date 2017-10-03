pins=7;
theta0=10;
chro=100;
simnum=10000000;
u=2.5*10^(-4);
msdata=zeros(simnum,pins);
%Ym=random('logn',N,N,[1,simnum]);
%Y=4.*Ym.*u;
Y=100.*rand(1,simnum);
for i=1:simnum 
    argv=sprintf('./ms %d 1 -t %f >msout',chro,Y(i));
    system(argv);
    infile=fopen('msout');
   % A=fscanf(infile,'%d',3);
   % B=fscanf(infile,'%c',5);
    fgets(infile);
    fgets(infile);
    fgets(infile);
    fgets(infile);
    num=fscanf(infile,'segsites:%d');
    if num==0
        Y(i)=0;
        msdata(i,:)=0;
	fclose(infile);
        continue;
    end
    data=zeros(chro,num);
  %  fscanf(infile,'%C',1);
    fgets(infile);
    fgets(infile);
    for ii=1:chro
        data(ii,:)=fscanf(infile,'%1d',num);
    end
    sfs=sum(data,1);
    for iii=1:num
        if sfs(iii)>num/2
            sfs(iii)=chro-sfs(iii);
        end
    end
    pin1=size(sfs(sfs<8),2);
    pin2=size(sfs(sfs<16),2)-pin1;
    pin3=size(sfs(sfs<24),2)-pin1-pin2;
    pin4=size(sfs(sfs<32),2)-pin1-pin2-pin3;
    pin5=size(sfs(sfs<40),2)-pin1-pin2-pin3-pin4;
    pin6=size(sfs(sfs<48),2)-pin1-pin2-pin3-pin4-pin5;
    pin7=num-pin1-pin2-pin3-pin4-pin5-pin6; 
    msdata(i,:)=[pin1 pin2 pin3 pin4 pin5 pin6 pin7];
    fclose(infile);
end
data=sum(msdata)/simnum;
save data.mat data;
save msdata.mat msdata;
save Y.mat Y;
