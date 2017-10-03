clear all;
pins=7;

chro=100;
simnum=500000;
u=2.5*10^(-4);
tdata=zeros(simnum,pins);  %%It is all called msdata, but in differnt model it may contain SFS or SFS+HFS summary statistics

Y=20.*rand(1,simnum);
parfor i=1:simnum 
    outfile=['msout',num2str(i)];
    argv=sprintf('./ms %d 1 -t %f >%s',chro,Y(i),outfile);
    system(argv);
    infile=fopen(outfile);
   % A=fscanf(infile,'%d',3);
   % B=fscanf(infile,'%c',5);
    fgets(infile);
    fgets(infile);
    fgets(infile);
    fgets(infile);
    num=fscanf(infile,'segsites:%d');
    if num==0
        Y(i)=0;
        tdata(i,:)=0;
        fclose(infile);     
        argv=sprintf('rm %s',outfile);
        system(argv);
        continue;
    end
    data=zeros(chro,num);
  %  fscanf(infile,'%C',1);
    fgets(infile);
    fgets(infile);
    for ii=1:chro
        data(ii,:)=fscanf(infile,'%1d',num);
    end
    fclose(infile);
    sfs=sum(data,1);
    pin1=size(sfs(sfs<8),2);
    pin2=size(sfs(sfs<16),2)-pin1;
    pin3=size(sfs(sfs<24),2)-pin1-pin2;
    pin4=size(sfs(sfs<32),2)-pin1-pin2-pin3;
    pin5=size(sfs(sfs<40),2)-pin1-pin2-pin3-pin4;
    pin6=size(sfs(sfs<48),2)-pin1-pin2-pin3-pin4-pin5;
    pin7=num-pin1-pin2-pin3-pin4-pin5-pin6; 
    tdata(i,:)=[pin1 pin2 pin3 pin4 pin5 pin6 pin7];

    argv=sprintf('rm %s',outfile);
    system(argv);
end
tdataAVE=sum(tdata)/simnum;
tY=Y';
save tdata.mat tdata;
save tY.mat tY;
