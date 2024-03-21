% Plot W-band reflectivity on 2009 November 26
% when the inversion jumped up and got thicker in the soundings
% after plot_1min_stat.m
read_parameters;
load([way_proc_data_wband '1min_stat/Z_1min.mat']);

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

%synchronize ceilo cloud base height clb10 to cloudbase10(time10yday)
[ii,loc]=ismember(round(clb10.yday*60*24/10),round(time10yday*60*24/10));
cloudbase10=NaN+time10yday;
cloudbase10(loc(ii))=clb10.height(ii);
% also sync ceilometer cloud fraction ceilfrac10(time10yday)
ceilfrac10=NaN+time10yday;
ceilfrac10(loc(ii))=clb10.cloudfrac(ii);

cloudbase10(cloudbase10>cloudtop10)=NaN;
cloudtop10(cloudbase10>cloudtop10)=NaN;

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

% load sondes
load([way_raw_data_balloon(1:end-4) 'mat/wtec_edt_interp2z.mat'])
yday_sonde=datenum(year,month,day,hour,minute,0)-datenum('0-jan-2008');

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
%colormap(drizzle) % problem with hot is can't discern -30dBZ from the blacked-out background & noise.

% PLOT timeheight of reflectivity and its variance
lines=false;
iiw=time_yday>=331 & time_yday<=332;
set(gcf,'color','k','inverthardcopy','off')
clf
ax(1)=subplot(2,1,2,'align');
% sondes
iis=182:187;
if lines; plot(yday_sonde([1;1;1;1;1]*iis),[-0.07;0;NaN;2;2.07]*(1+0*iis),'y-','clipping','off'); end;
hold on
b2rcolormap(21);
pcolor(time_yday(iiw),(Z.height-12.5)/1e3,Z.plotmean(iiw,:)'); shading flat; % moved down 1/2 range gate
plot(time10yday,lwp10/1e3,'c.','markersize',2.5) % LWP kg/m^3
set(ax(1),'ydir','normal','tickdir','out','ylim',[0 2],'color',[0 0 0],'xtick',316:2:334,'xminortick','on','fontsize',14)
set(ax(1),'xcolor','w','ycolor','w')
axis([331 332 0 2])
set(gca,'xtick',331:1/12:332,'xticklabel',0:2:24)
caxis([-30 10])
box on
title('Mean reflectivity (dBZ)','color','w')
ylabel('height (km)')
cb(1)=colorbar; set(cb(1),'fontsize',14)
% cloud base and top
plot(ax(1),time10yday,cloudtop10/1e3,'g.','markersize',2.5)
plot(ax(1),time10yday,cloudbase10/1e3,'m.','markersize',2.5)
xlabel('November 26 hour')

ax(2)=subplot(2,1,1,'align');
pcolor(yday_sonde-2/24,hz(hz>=800 & hz<=2100)/1e3,tz(hz>=800 & hz<=2100,:)); shading flat;
hold on
set(ax(2),'xcolor','w','ycolor','w','fontsize',14)
axis([329 334 .8 2.1])
set(gca,'tickdir','out','xtick',[329:330 331:1/12:332 333:334],'xticklabel',{'24' '25' '26','','','','','','','','','','','','27' '28' '29'})
cb(2)=colorbar; set(cb(2),'fontsize',14)
if lines; plot(yday_sonde([1;1;1;1;1]*iis),[0.77;0.8;NaN;2.1;2.13]*(1+0*iis),'y-','clipping','off'); end;
title('Temperature (C)','color','w')
ylabel('height (km)')
xlabel('2008 November date')
if lines
    print('-dpng',[way_proc_images_wband 'sond_refl_nov26.lines.png'])
    print('-depsc2',[way_proc_images_wband 'sond_refl_nov26.lines.eps'])
else
    print('-dpng',[way_proc_images_wband 'sond_refl_nov26.png'])
    print('-depsc2',[way_proc_images_wband 'sond_refl_nov26.eps'])
end
