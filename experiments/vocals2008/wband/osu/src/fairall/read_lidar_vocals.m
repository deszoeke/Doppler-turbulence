disp('read lidar vocals')
inDir='~/Data/cruises/VOCALS_2008/RHB/lidar/Processed/';

startDateNum = ddd+ehr/24+datenum([2007,12,31,0,0,0]);%ddd=julian day, ehr=hour
stopDateNum = ddd+(ehr+1)/24+datenum([2007,12,31,0,0,0]);
[out_zProf] = readZProfProfilesTxt(inDir,startDateNum,stopDateNum);
[out_profData] = readSDIRProfilesTxt(inDir,startDateNum,stopDateNum);
[n m]=size(out_zProf);
L=250;length(out_zProf(1).height)-10;
htl=[];jdl=[];LSNR=[];wvar=[];uspd=[];sgspd=[];x=[];j=1;
for i=1:m;
    m1=length(out_zProf(i).height);
    if m1>=L
        xx=out_zProf(i).decTime;
        jdl(j)=xx-datenum([2007,12,31,0,0,0]);
        x=out_zProf(i).height;
        htl(j,1:L)=x(1:L);
        x=out_zProf(i).mnSNR;
        LSNR(j,1:L)=x(1:L);
        x=out_zProf(i).atmosVarNsSpec;
        wvar(j,1:L)=x(1:L);
        x=out_profData(i).z_alt;
        htu(j,1:L)=x(1:L);
        x=out_profData(i).windSpeed;
        uspd(j,1:L)=x(1:L);
        x=out_profData(i).windStd;
        sgspd(j,1:L)=x(1:L);
        j=j+1;
    end;
end;
m=j-1;
for i=1:m
        ii=find(diff(htl(i,:))>40);
        if ~isempty(ii)
            hx(i)=htl(i,ii);
            hy(i)=htl(i,ii+1);
        else
            hx(i)=NaN;
            hy(i)=NaN;
        end;
end;

htz=30:30:2000;
LSNRm=[];wvarm=[];uspdm=[];sgspdm=[];

for i=1:m
    wvarm(i,:)=interp1(htl(i,:),wvar(i,:),htz);
    LSNRm(i,:)=interp1(htl(i,:),LSNR(i,:),htz);
    ii=find(LSNRm(i,:)==max(LSNRm(i,:)));
    wvarm(i,ii-1:66)=NaN;
    ij=find(wvarm(i,:)>5 | wvarm(i,:)<0 | LSNRm(i,:)<0);wvarm(i,ij)=NaN;
end;

wvarmn=nanmean(wvarm);
figure;plot(wvarmn,htz);xlabel('\sigma _W^2 (m^2/s^2)');ylabel('Altitude (m)');
    

