
h=(cldtop(:,2));%cloudtop
kk=find(ht<h);
topi=max(kk);
cldthk=(cldtop(:,8))/12600*25/1e3;
kk=find(ht<h-cldthk);
kk=find(detectm>.2);%return bottom (20%)
boti=kk(1);
dkcld=min(30,topi-boti+2);%

%%%  linear regression dBZ vs fall velocity
j=1;
dbb=-35:0;
jj=find(Ref(topi-2,:)<-20);
clear Pi ygg;
for dk=3:dkcld;
    jj=find(Ref(topi-dk,:)>-35);
    [P,S] = polyfit(Ref(topi-dk,jj),DV(topi-dk,jj),2);
    yg=polyval(P,dbb);
    ygg(j,:)=yg;
    Pi(j,:)=P;
    j=j+1;
end
figure;plot(dbb,mean(ygg),dbb,polyval(mean(Pi),dbb),'o');%mean reqressio fit

j=1;
velk=[];
for dk=1:dkcld;
vel = polyval(mean(Pi),Ref(topi-dk,:));
velk(j,:)=vel;
j=j+1;
end

vc=median(DV(topi-dkcld+1:topi,:)-velk);

%Regression on exponential fit
dbt=[ones(size(dbb')) exp(.06*dbb')];
%Regression on mean dBZ vs W
tg=nanmean(Ref(topi-dkcld:topi-floor(dkcld/2),:))';   % average over ?lower half of cloud?
vg=nanmean(DV(topi-dkcld:topi-floor(dkcld/2),:))'-vc';
[Pg,Sg] = polyfit(tg,vg,2);
xp=[ones(size(tg)) exp(.06*tg)];
aa=xp\vg;
ygb=dbt*aa;
%figure; dock; clf
plot(tg,vg,'.','color',[0.75 0 0.75],'markersize',1); % empirical data
hold on
plot(dbb,mean(ygg),'-ob')                             % avg quadratic fit Pi
plot(dbb,polyval(Pg,dbb),'--r')                       % lower half of cloud quadratic fit Pg
plot(dbb,ygb,'-x','color',[0 0.75 .75])               % exponential fit
xlabel('dBZ');ylabel('Vertical Motion (m/s)');

velk2=[];velk3=[];
j=1;
for dk=1:dkcld;
    vel = polyval(Pg,Ref(topi-dk,:));
    velk2(j,:)=vel;
    tgg=Ref(topi-dk,:)';xpg=[ones(size(tgg)) exp(.06*tgg)];
    velg=xpg*aa;
    velk3(j,:)=velg';
    j=j+1;
end

dk=5;vc2=mean(DV(topi-dkcld:topi-dk,:)-velk2(dk:dkcld,:));vc3=mean(DV(topi-dkcld:topi-dk,:)-velk3(dk:dkcld,:));

wp=[.02];ws= [.01];rp=3;rs=15;[n,wmn]=cheb1ord(wp,ws,rp,rs);
[b,a]=cheby1(n,rp,wmn);
vcfL=filter(b,a,vc-mean(vc));%\
vcfH=vc-vcfL-mean(vc);
vcfL3=filter(b,a,vc3-mean(vc3));%\
vcfH3=vc3-vcfL-mean(vc3);

if ~exist(DV_mot) % motion should be loaded already from WbandPlot
    % read_Wband_motion
    % DV_mot=interp1(hhk,vkons(9,:),hh);
    % ii=find(isnan(DV_mot));DV_mot(ii)=0;
    
    % Read Kongsberg motion compensation file
    kongfile=dir([way_raw_data_wband 'motion_adjT/' sprintf('%04i%03i%02i',[year,ddd,jam]) '*Kongsberg_adjT.txt']);
    kongfilename=[way_raw_data_wband 'motion_adjT/' kongfile(1).name];
    [kongtime,kongw,kongpitch,kongroll]=read_kongsberg(kongfilename);
    % synchronize Kongsberg and radar time (yearday)
    base_time_mld=double(base_time)/86400 + datenum(1970,1,1,0,0,0);
    base_time_yday=base_time_mld-datenum(year,1,0,0,0,0);
    [junk,ii,jj]=unique(kongtime);
    DV_mot=interp1(kongtime(ii),kongw(ii),base_time_yday+time_offset/86400)';
end

figure;
km=floor(dkcld/2);
ii=find(Ref(topi-km,:)>-50);
subplot(3,1,1);plot(hh(ii),DV(topi-km,ii)-velk2(km,ii),hh(ii),DV(topi-km,ii));axis([hh(100) hh(1200) -2 2]);ylabel('Wvelocity (m/s)');
subplot(3,1,2);plot(hh(ii),DV_mot(ii));axis([hh(100) hh(1200) -2 2]);ylabel('Wmotion (m/s)');
subplot(3,1,3);plot(hh(ii),DV(topi-km,ii)-velk2(km,ii)+DV_mot(ii));hold on;
plot(hh(ii),DV(topi-km-1,ii)-velk2(km-1,ii)+DV_mot(ii),'r');
plot(hh(ii),DV(topi-km-2,ii)-velk2(km-2,ii)+DV_mot(ii),'g');
axis([hh(100) hh(1200) -2 2]);ylabel('Wvelocity (m/s)');xlabel('Time (hrs)');
hold off;

fw=1/.3;

[S_w,F]=psd2(detrend(DVa(topi-km,:)),length(hh),fw);
[S_ws,Fsl]=specsmoo(S_w,fw);
[S_w2,F]=psd2(detrend(DV(topi-km,ii)-velk2(km,ii)),length(hh),fw);
[S_ws2,Fsl]=specsmoo(S_w2,fw);
lss=length(Fsl);
dff=diff(Fsl);
dff=[dff(1) dff];
Snoise=min(S_ws(lss-10:lss));
af=mean(Fsl(lss-38:lss-30).^(5/3).*(S_ws2(lss-38:lss-30)-Snoise));
figure;loglog(Fsl,(S_ws-Snoise).*Fsl,Fsl,S_ws2.*Fsl,Fsl(lss-41:lss),af*Fsl(lss-41:lss).^(-2/3));
title([ 'VOCALS Wband Radar Spectrum:   Day ',num2str(ddd),'   Hour ',num2str(hh(1))]);
xlabel('Frequency in Hz');
ylabel('W-Spectrum * Frequency');

read_lidar_vocals;

figure;plot(hh,vcfL);
kk=find(htu(1,:)>100);ux=uspd(1,kk(1));%lidar wind speed at 100m
j=1;Swi=[];kss=25;eps=[];Swif=[];sigvel=[];skewvel=[];hep=[];velg=[];Snoisez=[];sigvelsp=[];
dkcld2=50;
for dk=1:dkcld2;
    hep(j)=ht(topi-dk);%height
    retfrc(j)=detectm(topi-dk);%fraction of time there is return
    %veldbz = polyval(Pg,Ref(topi-dk,:));%estimate particle fall velocity
    tgg=Ref(topi-dk,:)';xpg=[ones(size(tgg)) exp(.06*tgg)];
    velq=xpg*aa;
    veldbz=velq';
    velf=DV(topi-dk,:)-veldbz;%Calc air motion
    velg=velf;
    jj=find(Ref(topi-dk,:)+39-20*log10(ht(topi-dk)')<8);%check for good SNR
    velf(jj)=0;velg(jj)=NaN;
    sigvel(j)=nanstd(velg);%sigmaW of good data
    skewvel(j)=nanmean((velg-nanmean(velg)).^3)/sigvel(j)^3;%skeweness W
    [S_w,F]=psd2(detrend(velf),length(hh),fw);%power spectrum
    [S_ws,Fsl]=specsmoo(S_w,fw);%Smoothed spectrum
    Snoise=min(S_ws(lss-10:lss));%estimate white noise level from hi freq region
    Snoisez(j)=Snoise;
    Swi(j,:)=S_ws;
    Swif(j,:)=(S_ws-Snoise).*Fsl;
    vls=.75/.54*(2*pi/ux)^.667*mean(Fsl(lss-38:lss-34).^(5/3).*(S_ws(lss-38:lss-34)-Snoise)); %used 6 points to estimate dissiaption
    eps(j)=vls^1.5;
    sigv2=sum((S_ws-Snoise).*dff);
    sigvelsp(j)=sigv2;
    j=j+1;
end
kk=find(retfrc>.5);%theshold on fraction of good signal
figure;
subplot(1,2,1);semilogx(eps(kk),hep(kk));title([ 'TKE Dissipation Rate']);axis([ 3.5e-5 3.5e-3 0 1.5]);
xlabel('\epsilon (m^2/s^3)');
ylabel('ALtitude (km)');
subplot(1,2,2);plot(sigvel(kk).^2,hep(kk),sigvelsp(kk),hep(kk),wvarmn,htz/1e3);title([ 'Variance:   Day ',num2str(ddd),'   Hour ',num2str(hh(1))]);axis([ 0 1 0 1.5]);
xlabel('\sigma ^2 (m^2/s^2)');
ylabel('ALtitude (km)');

% figure;loglog(Fsl,Swif);axis([1e-4 10 1e-4 1]);title([ 'VOCALS Wband Radar Spectrum:   Day ',num2str(ddd),'   Hour ',num2str(hh(1))]);
% xlabel('Frequency in Hz');
% ylabel('W-Spectrum * Frequency');
% 

