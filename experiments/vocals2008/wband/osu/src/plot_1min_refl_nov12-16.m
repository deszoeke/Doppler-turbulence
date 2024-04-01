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

ncpath='/Users/sdeszoek/Data/cruises/VOCALS_2008/RHB/radar/wband/Processed/cloudheight/';
% note netcdflib does not accept ~/ as a path
ncf=[ncpath 'VOCALS2008CloudBoundaries10min_0_2.nc'];
yday=nc_varget(ncf,'yday');
lat=nc_varget(ncf,'lat');
lon=nc_varget(ncf,'lon');
cloudtop=nc_varget(ncf,'cloudtop');
sondecloudtop=nc_varget(ncf,'sondecloudtop');
sondeinterpflag=nc_varget(ncf,'sondeinterpflag');
cloudbase=nc_varget(ncf,'cloudbase');
ncloud=nc_varget(ncf,'numceilcloud');
time10yday=yday;

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

% load sondes
load([way_raw_data_balloon(1:end-4) 'mat/wtec_edt_interp2z.mat'])
yday_sonde=datenum(year,month,day,hour,minute,0)-datenum('0-jan-2008');

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
max(Z.max,[],2)

% Column maximum reflectivity (CMR)
% vertical max reflectivity of 1-minute means
[zmx,k]=max(Z.mean,[],2);
% make 10-min max dBZ
zmx10=NaN+time10yday;
for ti=1:length(time10yday)
    ii=Z.time_yday>=time10yday(ti) & Z.time_yday<time10yday(ti)+10/60/24;
    if sum(ii)>0
        zmxmx10(ti)=max(zmx(ii));
        zmxmd10(ti)=nanmedian(zmx(ii));
    end
end

% histogram of 1-minute vertical max dBZ
edges=-40:39;
count=histc(zmx,edges);
count10=histc(zmxmd10,edges);
lwpbin=NaN+edges;
lwpbinx=lwpbin;
lwpbinn=lwpbin;
for bin=1:length(edges)-1;
    ii=zmxmd10>=edges(bin) & zmxmd10<edges(bin+1);
    if sum(ii)>0
        lwpbin(bin)=nanmedian(lwp10(ii));
        lwpbinx(bin)=nanmax(lwp10(ii));
        lwpbinn(bin)=nanmin(lwp10(ii));
    end
end

ax(3)=subplot(3,1,3);
hbar=bar(edges,count/60,'histc');
hold on
plot(edges+0.5,lwpbin/10,'k.');
plot(edges+0.5,lwpbinx/10,'k-');
plot(edges+0.5,lwpbinn/10,'k-');

set(hbar,'facecolor',0.7*[1 1 1])
axis([-40 25 0 55])
set(gca,'ytick',0:10:50)
set(gca,'tickdir','out','fontsize',14)
xlabel('max dBZ')
ylabel('frequency (hours)')
xlabel('column max reflectivity (dBZ)') 
set(gca,'tickdir','out')

% PLOT timeheight of reflectivity and its variance
iiw=time_yday>=317 & time_yday<321; % Nov 12-15.9999
iis=yday_sonde>=317 & yday_sonde<321;
ax(1)=subplot(3,1,1,'align');
pcolor(time_yday(iiw),(Z.height-12.5)/1e3,Z.plotmean(iiw,:)'); shading flat; % moved down 1/2 range gate
hold on
caxis([-35 15])
colormap(1-gray(20));
set(ax(1),'tickdir','out','fontsize',14)
axis([317 321 0 1.6])
set(gca,'xtick',317:321,'xticklabel',12:16)
title('Mean reflectivity (dBZ)')
ylabel('height (km)')
%cb(1)=colorbar; set(cb(1),'fontsize',14)
% cloud base and top
plot(ax(1),time10yday,cloudtop10/1e3,'k')
plot(ax(1),time10yday,cloudbase10/1e3,'k')
%plot(C.time_yday,C.prain/1e3)
%plot(Z.time_yday,max(Z.max,[],2)/20+.5,'m')

ax(2)=subplot(3,1,2,'align');
plot(time10yday,lwp10,'k') % LWP kg/m^3
hold on
plot(time10yday,adlwp10,'color',0.7*[1 1 1],'markersize',2.5) % LWP kg/m^3
%plot([317 321],[0 0],'k:')
set(ax(2),'fontsize',14)
axis([317 321 0 320])
set(gca,'xtick',317:321,'xticklabel',12:16,'tickdir','out')
ylabel('LWP (g m^{-2})')
xlabel('2008 November date')
    caxis([-35 15])
    colormap(1-gray(20));
    cb(1)=colorbar('northoutside'); set(cb(1),'fontsize',14)
% was called sond_refl_nov12-16.png, sond_refl_nov12-16.eps
% print('-dpng',[way_proc_images_wband 'refl_lwp_nov12-16.png'])
% print('-depsc2','-painters',[way_proc_images_wband 'refl_lwp_nov12-16.eps'])

% new adiabatic LWP
plot(time10yday,0.001*(cloudtop-cloudbase).^2)
% deviation from adiabatic, infer due to drizzle
plot(time10yday,lwp10-0.001*(cloudtop-cloudbase).^2,'color',[0 .7 0])
% predict LWP from max dBZ regression
plot(time10yday,300+9.2*zmxmd10,'r')
orient tall

%print('-dpng',[way_proc_images_wband 'refl_lwp_nov12-16_4PZ.png'])
%print('-depsc2','-painters',[way_proc_images_wband 'refl_lwp_nov12-16_4PZ.eps'])
