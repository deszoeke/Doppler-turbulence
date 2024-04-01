% Plot 10-minute statistics of Z and w from W-band radar.
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

% W-band cloud fraction from different thresholds
% see proc_wband_cloudfraction.m, proc_wband_cloudsfrac_10min.m
ncf2=[way_proc_data_wband 'cloudfraction/VOCALS2008_WbandCloudFraction10min_1_0a.nc']; % note netcdflib does not accept ~/ as a path
wcf.yday=nc_varget(ncf2,'yday'); % identical to cth10.yday
wcf.thr=nc_varget(ncf2,'reflectivity_threshold');
wcf.cft=nc_varget(ncf2,'cloud_fraction_threshold');

% load 10-min CL31 cloud base height and cloud fraction
A=load('/Users/sdeszoek/Data/cruises/VOCALS_2008/RHB/ceilometer/Processed/VOCALS2008ceilo_10min_a.txt');
clb10.yday=A(:,1);
clb10.cloudfrac=A(:,9);
clb10.height=A(:,13); % 85 percentile 1st cloud base height
clb10.height(clb10.height>2000)=NaN; % filter above 2000m
% load cb, ceilometer-bacscatter threshold cloud fraction (10-min)
load([way_proc_data_ceilo 'cloudfrac_backscatterthreshold.mat']);
%cb.yday == wcf.yday

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

% load inversion heights and superinversion thicknesses
invfile=[way_proc_data_balloon 'VOCALS2008_sonde_inversion.nc'];
sonde_yday=nc_varget(invfile,'yday');
sonde_lat=nc_varget(invfile,'lat');
sonde_lon=nc_varget(invfile,'lon');
hinvtop=nc_varget(invfile,'height_inv_top');
hinvbase=nc_varget(invfile,'height_inv_base');
tinvtop=nc_varget(invfile,'temp_inv_top');
tinvbase=nc_varget(invfile,'temp_inv_base');
thinvtop=nc_varget(invfile,'theta_inv_top');
thinvbase=nc_varget(invfile,'theta_inv_base');
dthdzinv=(thinvtop-thinvbase)./(hinvtop-hinvbase);

% synchronize times (1 min)
temp=time_yday(1)+0:1/(60*24):time_yday(end);  % even 1-min time_yday without gaps
time1yday=NaN+zeros(10*ceil(length(temp)/10),1); % pad w/ nans to get length multiple of 10
time1yday(1:length(temp))=temp;
cloudfrac1=NaN+zeros(length(time1yday),1);
[ii,loc]=ismember(round(time_yday*60*24),round(time1yday*60*24));
% time1yday(loc(ii)) approx. equals time_yday(ii)
cloudfrac1(loc(ii))=cloud_fraction(ii);
% 10-minute median of 1-minute data
time10yday=time1yday(1:10:end); % even 10-min
n10m=length(time10yday);
cloudfrac10a=nanmean(reshape(cloudfrac1,[10 n10m]))';

% synchronize LWP and Ztopmax10(time10yday)
[ii,loc]=ismember(round(lwp.yday*60*24/10),round(time10yday*60*24/10));
% time10yday(loc(ii)) approx equals lwp.yday(ii)
lwp10=NaN+time10yday;
adlwp10=lwp10;
lwp10(loc(ii))=lwp.lwp(ii);
adlwp10(loc(ii))=lwp.adlwp(ii);

% synchronize W-band cloud top height cth10 to cloudtop10(time10yday)
[ii,loc]=ismember(round(cth10.yday*60*24/10),round(time10yday*60*24/10));
cloudtop10=NaN+time10yday;
cloudtop10(loc(ii))=cth10.height(ii);
cloudfrac10=NaN+time10yday;
cloudfrac10(loc(ii))=cth10.frac(ii);
cloudfracth10=NaN+zeros(max(loc),length(wcf.thr));
cloudfracth10(loc(ii),:)=wcf.cft(ii,:);

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
Z.plotmax=Z.max;
Z.plotmax(isnoise)=NaN;
w.plotmean=w.mean;
w.plotmean(isnoise)=NaN;

plot(time_yday,nanmax(Z.plotmax,[],2))
hold on
plot(sonde_yday,1e3*dthdzinv,'color',[0 .5 0])

iloc=sonde_lat<-15 & sonde_lon<-72;
hist(dthdzinv(iloc)*1e3,40); % max is ~260, K/km wide
set(gca,'fontsize',16)
xlabel('d\theta/dz (K/km)')
ylabel('frequency')
ylabel('number of sondes')
set(gca,'tickdir','out','fontsize',16)
print('-dpng',[way_proc_images_balloon 'inv_stability.png'])
