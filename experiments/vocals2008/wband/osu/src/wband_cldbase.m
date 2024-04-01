read_parameters;

% load 1-minute cloud fraction and top height data (Simon de Szoeke, Oregon State U.)
base_time=1225916508;
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
time_yday=time_top_yday;

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

% synchronize LWP and time10yday
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


% load height-resolved 1-min Wband reflectivity Z
load([way_proc_data_wband '1min_stat/Z_1min.mat']);
% cloud base from max reflectivity (along vertical)
[zmx,k]=max(Z.mean,[],2);
% height of max reflectivity
hzmx=Z.height(k);
hzmx(zmx<-20)=NaN; % eyeball test of whether cloud is precipitating

pcolor(Z.time_yday,Z.height,max(Z.mean',-40)); shading flat;
caxis([-40 20])
hold on
plot(Z.time_yday,hzmx,'m-')
colorbar
% ...compare to ceilometer

edges=-40:39;
count=histc(zmx,edges);
bar(edges,count,'histc')
axis([-40 40 0 3200])
set(gca,'fontsize',18)
title('1-min max dBZ histogram')
ylabel('count')
xlabel('max dBZ') 
print('-depsc2',[way_proc_images_wband 'zmax_hist.eps'])

% average 600 m dBZ
hi=find(Z.height>600,1,'first'); % 600 m height index
for ti=1:length(time10yday)
    ii=Z.time_yday>=time10yday(ti) & Z.time_yday<time10yday(ti)+10/60/24;
    temp=Z.mean(ii,hi+[-1 0 1]);
    zmxmed=nanmedian(zmx(ii));
    dBZ600(ti)=nanmean(temp(:));
    dBZmed600(ti)=nanmedian(temp(:));
    dBZ85600(ti)=quantile(temp(:),0.85);
    fracgtm10(ti)=sum(temp(:)>-10)/sum(isfinite(temp(:)));
    fracgtm5(ti)=sum(temp(:)>-5)/sum(isfinite(temp(:)));
    fracgt0(ti)=sum(temp(:)>0)/sum(isfinite(temp(:)));
    fracgt5(ti)=sum(temp(:)>5)/sum(isfinite(temp(:)));
    fracgt10(ti)=sum(temp(:)>10)/sum(isfinite(temp(:)));
end
% 85th percentile dBZ is off from 50th by a nearly constant offset ~5dB
% almost no 10 minute samples have dBZ>0
% 1.96% are >-10 dBZ
% 1.14% are >-5 dBZ
% 0.57% are >  0 dBZ
% 0.32% are >  5 dBZ
% 0.13% are > 10 dBZ