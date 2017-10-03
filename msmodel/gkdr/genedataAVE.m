%function [ obsdata] = genedata(theta)
clear all;
pins=7;
simnum=5000;
chro=100;
obsdata=zeros(simnum,pins);  %%It is all called msdata, but in differnt model it may contain SFS or SFS+HFS summary statistics
%Ym=random('logn',N,N,[1,simnum]);
%Y=4.*Ym.*u;
Y=15;
parfor i=1:simnum


    outfile=['msout',num2str(i)];
    argv=sprintf('./ms %d 1 -t %f >%s',chro,Y,outfile);
    system(argv);
    infile=fopen(outfile);
   % A=fscanf(infile,'%d',3);
   % B=fscanf(infile,'%c',5);
    fgets(infile);
    fgets(infile);
    fgets(infile);
    fgets(infile);
    num=fscanf(infile,'segsites:%d');
    temdata=zeros(chro,num);

  %  fscanf(infile,'%C',1);
    fgets(infile);
    fgets(infile);
    for ii=1:chro
        temdata(ii,:)=fscanf(infile,'%1d',num);
    end
    fclose(infile);
    sfs=sum(temdata,1);
    pin1=size(sfs(sfs<8),2);
    pin2=size(sfs(sfs<16),2)-pin1;
    pin3=size(sfs(sfs<24),2)-pin1-pin2;
    pin4=size(sfs(sfs<32),2)-pin1-pin2-pin3;
    pin5=size(sfs(sfs<40),2)-pin1-pin2-pin3-pin4;
    pin6=size(sfs(sfs<48),2)-pin1-pin2-pin3-pin4-pin5;
    pin7=num-pin1-pin2-pin3-pin4-pin5-pin6; 
    obsdata(i,:)=[pin1 pin2 pin3 pin4 pin5 pin6 pin7];
    
    argv=sprintf('rm %s',outfile);
    system(argv);
end
obsdataAve=sum(obsdata)/simnum;
save ObAVE.mat obsdataAve;
%end
