% Plot 1-minute statistics of Z and w from W-band radar.
% VOCALS 2008 :: 2009-07-07 :: Simon de Szoeke

read_parameters;

load([way_proc_data_wband '1min_stat/Z_1min.mat']);
load([way_proc_data_wband '1min_stat/w_1min.mat']);

% Define dimensions
% vertical coordinate (center of ranges, m)
% h=load([way_proc_data_wband '1min_stat/height.txt']);
% % dBZ and w bin edges
% Zbins=[-inf -40:2:30 inf];
% wbins=[-inf -5:0.2:5 inf];
% % time_aggreg in seconds since base_time
base_time=1225916508;
% time_yday=datenum(0,0,0,0,0,base_time+time_aggreg)+datenum(1970,1,1,0,0,0)-datenum(2008,1,0,0,0,0); %yearday

nheight=length(Z.height); % 75
nmin=60;
nhour=24;
nZbin=length(Z.bins);
nwbin=length(w.bins);

% load C-band 3-min tilt 1 25-60 km annulus reflectivity histograms (Sandra Yuter, North Carolina State U.)
A=load('~/Data/cruises/VOCALS_2008/RHB/radar/C-band/leg2nocoasts.cz.tilt1.25to60km.stats');
C.year=A(:,1);
C.month=A(:,2);
C.date=A(:,3);
C.hour=A(:,4);
C.minute=A(:,5);
C.minim=A(:,6);
C.maxim=A(:,7);
C.mean=A(:,8);
C.hist=A(:,9:end);
C.Zbin=-30:59;
C.time_yday=datenum(C.year,C.month,C.date,C.hour,C.minute,0)-datenum(2008,1,0,0,0,0);
C.prain=sum(C.hist(:,C.Zbin>10),2);

% load 1-minute cloud fraction and top height data (Simon de Szoeke, Oregon State U.)
A=load([way_proc_data_wband 'cloudheight/CloudHeight_1min_2008310-336.txt']);
A(A==-999)=NaN; % NaN out missing value of -999
time=A(:,1); % seconds since base_time
cloud_fraction=A(:,2);
cloud_top_mean=A(:,3);
cloud_top_std=A(:,4);
cloud_top_val15=A(:,5);
cloud_top_median=A(:,6);
cloud_top_val85=A(:,7);
cloud_top_valmean=A(:,8); %mean of tops within 15-85 percentile
num_returns=A(:,9);
num_cloud_returns=A(:,10);
num_clear_returns=A(:,11);
num_nonmet_returns=A(:,12);
flag=A(:,13);
% flag values
% 0 no errors detected, cloud top computed.
% 1 cloud fraction=0, cloud top set to NaN.
% 2 cloud fraction>=0 but found no 3 cloud returns with consistent height,
%   possible glitch, provisional cloud top computed.
% NaN out percentile-based cloud tops=0 using flag
lf=logical(flag);
cloud_top_val15(lf)=NaN;
cloud_top_val85(lf)=NaN;
cloud_top_valmean(lf)=NaN;
cloud_top_median(lf)=NaN;
time_top_yday=datenum(0,0,0,0,0,base_time+time)+datenum(1970,1,1,0,0,0)-datenum(2008,1,0,0,0,0);
% ==time_yday everywhere

% load 10-min W-band cloud top height (Simon de Szoeke, Oregon State U.)
chf='/Users/sdeszoek/Data/cruises/VOCALS_2008/RHB/radar/wband/Processed/cloudheight/VOCALS2008_WbandCloudHeight10min_1_0.nc';
cth10.yday=nc_varget(chf,'yday');
cth10.height=nc_varget(chf,'cloudtop');
cth10.frac=nc_varget(chf,'cloudfrac');

% load 10-min CL31 cloud base height and cloud fraction
A=load('/Users/sdeszoek/Data/cruises/VOCALS_2008/RHB/ceilometer/Processed/VOCALS2008ceilo_10min_a.txt');
clb10.yday=A(:,1);
clb10.cloudfrac=A(:,9);
clb10.height=A(:,13); % 85 percentile 1st cloud base height
clb10.height(clb10.height>2000)=NaN; % filter above 2000m

% load 10-min LWP (Paquita Zuidema, U. Miami)
lwpfile=[way_proc_data_umiami 'vocals_rhb_lwp10min_v1.nc'];
lwp.yday =nc_varget(lwpfile,'jd_jdfrac');
lwp.lwp  =nc_varget(lwpfile,'lwp'); % g/m^2, 2-band mailbox physical retrieval (Zuidema 2005)
lwp.nobs =nc_varget(lwpfile,'nlwp');
lwp.adlwp=nc_varget(lwpfile,'adlwp');
% QC data
lwp.lwp(lwp.nobs<10)=NaN;
lwp.adlwp(lwp.adlwp>1e10 | lwp.adlwp<0)=NaN;

% load flux file for position
A=load([way_proc_data_flux 'cat/flux_5hf_VOCALS2008.txt']);
flux.yday=A(:,1);
flux.lat=A(:,17);
flux.lon=A(:,18);

% line up all the cloud tops
[H,Y]=meshgrid(Z.height,cloud_top_val85);
ct_rel_height=H-Y;
[X,T]=meshgrid(Z.height,time_yday); % lines them up for a plot
if false % plot
    pcolor(T,ct_rel_height+nanmean(cloud_top_mean),Z.mean)
    shading flat
    set(gca,'xlim',[316 334.5],'ylim',[-400 2000],'color','k')
end

% get maximum reflectivity w/in 100 m of cloud top
% itop=abs(ct_rel_height)<10; % old way
% [ii,jj]=find(itop);
% plot(time_yday(ii),Z.height(jj),'.')
% Z_top_max=max(Z.mean(ii,jj-4:jj),[],2);
Zmask=Z.mean;
Zmask(ct_rel_height>10 & ct_rel_height<-110)=NaN;
test=max(Zmask,[],2);
test(test<-32)=NaN;
Ztopmax=10.^(test/10); % reflectivity Z from dBZ
% get cloud top vertical velocity
wmask=Zmask;

% synchronize times (1 min)
temp=time_yday(1)+0:1/(60*24):time_yday(end);  % even 1-min time_yday without gaps
time1yday=NaN+zeros(10*ceil(length(temp)/10),1); % pad w/ nans to get length multiple of 10
time1yday(1:length(temp))=temp;
Ztopmax1=NaN+zeros(length(time1yday),1);
cloudfrac1=NaN+zeros(length(time1yday),1);
[ii,loc]=ismember(round(time_yday*60*24),round(time1yday*60*24));
% time1yday(loc(ii)) approx. equals time_yday(ii)
Ztopmax1(loc(ii))=Ztopmax(ii);
cloudfrac1(loc(ii))=cloud_fraction(ii);
% 10-minute median of 1-minute data
time10yday=time1yday(1:10:end); % even 10-min
n10m=length(time10yday);
Ztopmax10=nanmedian(reshape(Ztopmax1,[10 n10m]))';
cloudfrac10a=nanmean(reshape(cloudfrac1,[10 n10m]))';
dBZtopmax10=10*log10(Ztopmax10);

% synchronize LWP and Ztopmax10(time10yday)
[ii,loc]=ismember(round(lwp.yday*60*24/10),round(time10yday*60*24/10));
% time10yday(loc(ii)) approx equals lwp.yday(ii)
lwp10=NaN+Ztopmax10;
adlwp10=lwp10;
lwp10(loc(ii))=lwp.lwp(ii);
adlwp10(loc(ii))=lwp.adlwp(ii);

% synchronize W-band cloud top height cth10 to cloudtop10(time10yday)
[ii,loc]=ismember(round(cth10.yday*60*24/10),round(time10yday*60*24/10));
cloudtop10=NaN+time10yday;
cloudtop10(loc(ii))=cth10.height(ii);
cloudfrac10=NaN+time10yday;
cloudfrac10(loc(ii))=cth10.frac(ii);

%synchronize ceilo cloud base height clb10 to cloudbase10(time10yday)
[ii,loc]=ismember(round(clb10.yday*60*24/10),round(time10yday*60*24/10));
cloudbase10=NaN+time10yday;
cloudbase10(loc(ii))=clb10.height(ii);
% also sync ceilometer cloud fraction ceilfrac10(time10yday)
ceilfrac10=NaN+time10yday;
ceilfrac10(loc(ii))=clb10.cloudfrac(ii);

cloudbase10(cloudbase10>cloudtop10)=NaN;
cloudtop10(cloudbase10>cloudtop10)=NaN;

% compute volume mixing ratio at top from LWP and cloud thickness
Qctop=2*lwp10./(cloudtop10-cloudbase10);
% mass mixing ratio
rhotop=1.0; % g m^-3, approx from rho = rhosfc exp(-1.4/8), rhosfc=1.19;
qctop=Qctop/rhotop;

% plot cloud top reflectivity vs mixing ratio
plot(Qctop,dBZtopmax10,'o','markersize',3)
set(gca,'xscale','log','fontsize',14)
hold on
plot(5e-1+[0 4.3/2],[-33 10],'r-')
axis([1e-1 1e1 -33 10])
title('cloud top reflectivity vs. mixing ratio');
ylabel('cloud top maximum dBZ');
xlabel('cloud volume mixing ratio Q_c (g m^{-3})');
%print('-dpng',[way_proc_images_wband 'maxz_Qc.png']);

% plot time series of dBZtopmax10 time10yday
ii=isfinite(time10yday+dBZtopmax10);
hl=area(time10yday(ii),dBZtopmax10(ii),-40);
set(hl,'edgecolor','none','facecolor',[1 .8 .7])
hold on
plot(time10yday,dBZtopmax10,'r.-','markersize',3,'linewidth',0.5)
axis([316 334.5 -34 15])
set(gca,'xminortick','on','xtick',316:2:334,'tickdir','out')
set(gca,'fontsize',14)
ylabel('cloud top maximum dBZ')
xlabel('2008 year day')
title('W-band cloud top maximum reflectivity')
%print('-dpng',[way_proc_images_wband 'zmax_time.png'])

% compute hourly median and mean diurnal cycles
local_UTC_offset=-80/360*24;
[cldtopmd,cldtopmn,cldtopn]=diurnal_hourly_median(local_UTC_offset+time10yday, cloudtop10);
[cldbasmd,cldbasmn,cldbasn]=diurnal_hourly_median(local_UTC_offset+time10yday, cloudbase10);
[cldthkmd,cldthkmn,cldthkn]=diurnal_hourly_median(local_UTC_offset+time10yday, cloudtop10-cloudbase10);
[dBZmaxmd,dBZmaxmn,dBZmaxn]=diurnal_hourly_median(local_UTC_offset+time10yday, dBZtopmax10);
[lwpmd   ,lwpmn   ,lwpn   ]=diurnal_hourly_median(local_UTC_offset+time10yday, lwp10);

% plot cloud diurnal cycles
plot(24*(mod(time10yday+local_UTC_offset,1)+rand(length(time10yday),1)/144),cloudtop10,'ko','markersize',3)
hold on
plot(24*(mod(time10yday+local_UTC_offset,1)+rand(length(time10yday),1)/144),cloudbase10,'r+','markersize',3)
plot(24*(mod(time10yday+local_UTC_offset,1)+rand(length(time10yday),1)/144),cloudtop10-cloudbase10,'b.','markersize',3)
plot(-0.5:1:24.5,cldtopmd([end 1:end 1]),'k','linewidth',1.3)
plot(-0.5:1:24.5,cldbasmn([end 1:end 1]),'r','linewidth',2)
plot(-0.5:1:24.5,cldthkmd([end 1:end 1]),'b','linewidth',1.3)
axis([0 24 0 2000])
set(gca,'xtick',0:6:24,'fontsize',14)
title('cloud height diurnal cycle')
xlabel('local hour')
ylabel('cloud top, base, and thickness (m)')
%print('-dpng',[way_proc_images_wband 'diurnal_cloud_height.png'])

% Cband rain proxy is subject to unusual cases, so does not lend itself to averaging or diurnal cycle comparison.
%plot(24*(mod(C.time_yday+local_UTC_offset,1)+rand(length(C.time_yday),1)/144),10*C.prain,'r.','markersize',2)
%hold on
%ha=area(-0.5:1:24.5,lwpmd([end 1:end 1])); set(ha,'facecolor',[.7 .9 1],'edgecolor','none')

% plot LWP diurnal cycle
plot(-0.5:1:24.5,lwpmd([end 1:end 1]),'linewidth',1.3);
hold on
plot(24*(mod(time10yday+local_UTC_offset,1)+rand(length(time10yday),1)/144),lwp10,'bo','markersize',3)
axis([0 24 0 600])
set(gca,'xtick',0:6:24,'fontsize',14)
title('liquid water path diurnal cycle')
xlabel('local hour')
ylabel('liquid water path (g m^{-2})')
%print('-dpng',[way_proc_images_wband 'diurnal_lwp.png'])

%diurnal wband max(Z)
plot(-0.5:1:24.5,dBZmaxmd([end 1:end 1]),'r-','linewidth',1.3);
hold on
plot(24*(mod(time10yday+local_UTC_offset,1)+rand(length(time10yday),1)/144),dBZtopmax10,'ro','markersize',3)
axis([0 24 -33 15])
set(gca,'xtick',0:6:24,'fontsize',14)
title('W-band cloud top max dBZ diurnal cycle')
xlabel('local hour')
ylabel('cloud top maximum reflectivity (dBZ)')
%print('-dpng',[way_proc_images_wband 'diurnal_maxz.png'])

% scatter LWP and Ztopmax10
plot(lwp10,dBZtopmax10,'.','markersize',3)
hold on
plot(3*[1 1000],[-50 10],'r') % power law max(Z)~LWP^2
set(gca,'yscale','linear','xscale','log','xlim',[1e1 8e2],'ylim',[-33 0])
set(gca,'xtick',[10 100],'xticklabel',[10 100],'fontsize',14)
ylabel('cloud-top max. reflectivity (dBZ)') % dBZ=10log_{10}[mm^6 m^{-3}/1 mm^6 m^{-3}]
xlabel('liquid water path (g m^{-3})')
title('Cloud-top max. reflectivity vs. LWP');
%print('-dpng',[way_proc_images_wband 'maxdBZ_LWP.png'])

% fit LWP and cloud geometry
ii=isfinite(cloudtop10-cloudbase10+lwp10);
P=polyfit(log(lwp10(ii)),log(cloudtop10(ii)-cloudbase10(ii)),1);
% scatter LWP and cloud geometry
clf
plot(lwp10,cloudtop10-cloudbase10,'b.','markersize',3)
hold on
plot(sort(lwp10),27*sort(lwp10).^0.5,'r-')
%plot(sort(lwp10),2*sort(lwp10),'r--')
axis([0 500 0 600])
set(gca,'xtick',0:100:500,'fontsize',14)
set(gca,'yscale','log','xscale','log')
title('cloud thickness vs. LWP')
xlabel('liquid water path (g m^{-3})')
ylabel('cloud thickness (m)');
%print('-dpng',[way_proc_images_wband 'cloudthick_LWP.png']);

% scatter cloud thickness vs. max(Z)
clf
plot(cloudtop10-cloudbase10,dBZtopmax10,'b.','markersize',3)
hold on
%h=sort(cloudtop10-cloudbase10);
%plot(h,10*log10(h.^2))
axis([0 600 -32 15])
set(gca,'xtick',0:100:500,'fontsize',14)
title('cloud top maximum reflectivity vs. cloud thickness')
ylabel('cloud top maximum dBZ')
xlabel('cloud thickness (m)');
%print('-dpng',[way_proc_images_wband 'cloudthick_zmax.png']);


% plot adiabatic vs. microwave LWP and cloud thickness
% (compares my new cloud thickness with what Paquita used)
plot(lwp10,adlwp10,'b^','markersize',3)
hold on
plot(cloudtop10-cloudbase10,adlwp10,'ro','markersize',3)
xlabel({'microwave LWP (g m^{-3})' 'cloud thickness (m)'})
ylabel 'adiabatic LWP (g m^{-3})'
axis equal
axis([0 600 0 400])
title('adiabatic LWP vs. microwave LWP (blue) and cloud thickness (red)')
%print('-dpng',[way_proc_images_wband 'LWP_adLWP_cloudthick.png']);

% plot cloud height versus max(Z)

% define a noise level for the reflectivity based on Ken Moran's manual
radar_const=19.613; % dB_
dB_offset=-1.72; % dB
min_detectable_signal=-113.304+dB_offset; % dBmW
analog_noise=-120.74;  % peak, dBmW
digital_noise=-115.4;   % peak, dBmW
threshold_margin=3.5; % dB
noisefloor=20*log10(Z.height)+radar_const+digital_noise+threshold_margin;
% do not plot colors/contours below the noise threshold
isnoise=Z.mean<repmat(noisefloor',[size(Z.mean,1) 1]);
Z.plotmean=Z.mean;
Z.plotmean(isnoise)=NaN;
Z.plotvar=Z.var;
Z.plotvar(isnoise)=NaN;
w.plotmean=w.mean;
w.plotmean(isnoise)=NaN;

% PLOT timeheight of reflectivity and its variance
set(gcf,'color','k','inverthardcopy','off')
ax(1)=subplot(2,1,1,'align');
colormap(hot(16)) % problem with hot is can't discern -30dBZ from the blacked-out background & noise.
pcolor(time_yday,(Z.height-12.5)/1e3,Z.plotmean'); shading flat; % moved down 1/2 range gate
hold on
plot(time10yday,lwp10/1e3,'c.','markersize',2.5) % LWP kg/m^3
set(ax(1),'ydir','normal','tickdir','out','ylim',[0 2],'color',[0 0 0],'xtick',316:2:334,'xminortick','on','fontsize',14)
set(ax(1),'xcolor','w','ycolor','w')
axis([316 334.5 0 2])
caxis([-30 10])
title('Mean reflectivity (dBZ)','color','w')
ylabel('height (km)')
cb=colorbar; set(cb,'fontsize',14)

ax(2)=subplot(2,1,2,'align');
pcolor(time_yday,(Z.height-12.5)/1e3,sqrt(Z.plotvar')); shading flat;
hold on
plot(C.time_yday,C.prain/1e3,'c.','markersize',2.5)
set(ax(2),'ydir','normal','tickdir','out','ylim',[0 2],'color',[0 0 0],'xtick',316:2:334,'xminortick','on','fontsize',14)
set(ax(2),'clim',[1 10],'xlim',[316 334.5],'ylim',[0 2])
set(ax(2),'xcolor','w','ycolor','w')
title('Standard deviation of reflectivity (dBZ)','color','w')
ylabel('height (km)')
xlabel('year day 2008')
cb=colorbar; set(cb,'fontsize',14)
%print('-dpng',[way_proc_images_wband 'wband_Z_leg2_lwp_prain.png'])
%print('-depsc2',[way_proc_images_wband 'wband_Z_leg2_lwp_prain.eps'])
plot(ax(1),time10yday,cloudtop10/1e3,'g.','markersize',2.5)
plot(ax(2),time10yday,cloudtop10/1e3,'g.','markersize',2.5)
%print('-dpng',[way_proc_images_wband 'wband_Z_leg2_lwp_prain_top.png'])
%print('-depsc2',[way_proc_images_wband 'wband_Z_leg2_lwp_prain_top.eps'])
plot(ax(1),time10yday,cloudbase10/1e3,'m.','markersize',2.5)
plot(ax(2),time10yday,cloudbase10/1e3,'m.','markersize',2.5)
%print('-dpng',[way_proc_images_wband 'wband_Z_leg2_lwp_prain_top_base.png'])
%print('-depsc2',[way_proc_images_wband 'wband_Z_leg2_lwp_prain_top_base.eps'])

% w looks noisy:
subplot(2,1,1,'align')
imagesc(time_yday,w.height/1e3,w.mean')
set(gca,'ydir','normal','tickdir','out','ylim',[0 2],'color',[0 0 0])
colormap([colormap(hot); 0 0 0])
title('Mean vertical velocity (m/s, +up)')
ylabel('height (km)')
colorbar
caxis([-1.5 1.5])
subplot(2,1,2,'align')
imagesc(time_yday,w.height/1e3,sqrt(w.var'))
set(gca,'ydir','normal','tickdir','out','ylim',[0 2],'color',[0 0 0])
title('Standard deviation of vertical velocity (m/s)')
xlabel('year day 2008')
colorbar
caxis([0 2])
%print('-dpng',[way_proc_images_wband 'wband_w_leg2.png'])

% cloud/drizzle colormap
drizzle= [ 0.3137    0.3176    0.3137
           0.4092    0.4118    0.4092
           0.5046    0.5059    0.5046
           0.6000    0.6000    0.6000
           0.6637    0.6608    0.6824
           0.7275    0.7216    0.7647
           0.7912    0.7824    0.8471
           0.8549    0.8431    0.9294
           0.6745    0.6059    0.9647
           0.4941    0.3686    1.0000
           0.2706    0.3608    1.0000
           0.2000    0.6000    1.0000
           0         1.0000         0
           1.0000    1.0000         0
           1.0000    0.8000         0
           1.0000    0.6000         0 ];
colormap(drizzle);

% look at November 26
set(gca,'xlim',[331 332])
set(gca,'ylim',[0 2])
set(gca,'xtick',331:1/12:332,'xticklabel',0:2:24)

% e.g. cfad
imagesc(Z.bins(2:end-1),Z.height,squeeze(Z.cfad(4000,:,2:end-1)));
set(gca,'ydir','normal')

% compare ceilometer and W-band cloud fraction, see plot_cloudfraction.m
plot(clb10.yday,clb10.cloudfrac,'.')
hold on
plot(time_top_yday,cloud_fraction,'r.')
axis([315 337 0 1])

clf
loglog(1-cloudfrac10,1-ceilfrac10,'r.')
xlabel('W-band')
ylabel('ceilometer')
title('clear fraction')
hold on
loglog(cloudfrac10,ceilfrac10,'.')
xlabel('W-band')
ylabel('ceilometer')
title('cloud fraction')

% cloud fraction histograms
cbin=[-Inf 0:0.01:1 Inf];
hist_ceilfrac=histc(ceilfrac10,cbin);
hist_cloudfrac=histc(cloudfrac10,cbin);
ceilfish=0.5*log((1+ceilfrac10)./(1-ceilfrac10));
cloudfish=0.5*log((1+cloudfrac10)./(1-cloudfrac10));

% TO DO: Effective radii or CCN implied by reflectivity, liquid water path, and cloud edges.
