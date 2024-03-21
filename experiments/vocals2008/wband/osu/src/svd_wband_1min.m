% svd_wband_1min.m THIS FILE MISNAMED!
% VOCALS 2008 :: 2010-06-19 :: Simon de Szoeke
%
% Singular value decomposition (empirical orthogonal function, EOF) decomposition of 
% W-band refelectivity (& w). Shift to align cloud top heights from 1-minute W-band radar.

% Requires data:
% cloudheight/CloudHeight_1min_2008310-336.txt from proc_wband_cloudtop.m
% 1min_stat/Z_1min.mat from compil_1min_stat.m

% preamble
% cd ~/Data/cruises/VOCALS_2008/RHB/Scientific_analysis/programs/VOCALS2008_programs/wband/
warning('off','MATLAB:interp1:NaNinY');
read_parameters

% load 1-min cloud fraction and top height data
A=load([way_proc_data_wband 'cloudheight/CloudHeight_1min_2008310-336.txt']);
A(A==-999)=NaN; % NaN out missing value of -999
time=A(:,1);
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

starter=1;
base_time(starter)=1225916508; % s since 1970-1-1 00:00
time_yday=datenum(0,0,0,0,0,base_time(starter)+time)+datenum(1970,1,1,0,0,0)-datenum(2008,1,0,0,0,0);

% load 1-min Z data
load([way_proc_data_wband '1min_stat/Z_1min.mat']);

% load 10-min cloud boundaries
cbfile=[way_proc_data_wband 'cloudheight/VOCALS2008CloudBoundaries10min_0_2.nc'];
cb_time=nc_varget(cbfile,'yday');
cb_height=nc_varget(cbfile,'cloudbase');

% load 10-min synthesis data (lat,lon, etc.)
xyz=load('~/Data/Stratus/cloud_10min2008.txt');
                    % variable name; units
yrr=xyz(:,1);       % Year of the cruise; Gregorian year
jd=xyz(:,2);        % yearday; days since proir December 31
lat=xyz(:,3);       % Latitude; degrees
lon=xyz(:,4);       % Longitude; degrees
tsg=xyz(:,5);       % Ship thermosalinograph water temperature (5 m depth); degrees C
ts=xyz(:,6);        % ESRL seasnake water temperature (5 cm depth); degrees Celsius
ta=xyz(:,7);        % Air Temperature (18 m); degree Celsius
ux=xyz(:,8);        % Westerly wind component; m/s
uy=xyz(:,9);        % Southerly wind component; m/s
u=xyz(:,10);        % Wind speed; m/s
hs=xyz(:,11);       % Sensible heat flux; W/m^2
hl=xyz(:,12);       % Latent heat flux; W/m^2
stress=xyz(:,13);   % Stress; N/m^2
rain=xyz(:,14);     % Rain rate; mm/hr
rl=xyz(:,15);       % Downward IR flux; W/m^2
rlclr=xyz(:,16);    % Clear sky downward IR flux; W/m^2
rs=xyz(:,17);       % Downward solar flux; W/m^2
rsclr=xyz(:,18);    % Clear sky downward solar flux; W/m^2
zb1m=xyz(:,19);     % 15th percentile cloud base height, 15% of cloud bases lower than zb1md; meters
zbm=xyz(:,20);      % Median cloud base height; meters
zb2m=xyz(:,21);     % 85th percentile cloud base height, 85% of cloud bases lower than zb2md; meters
top=xyz(:,22);      % Cloud top height; meters
clodthk=xyz(:,23);  % Cloud thickness; meters
tau=xyz(:,24);      % Cloud optical thickness; none
lwp=xyz(:,25);      % Cloud liquid water path; g/m^2
cf=xyz(:,26);       % Cloud fraction; none
nd=xyz(:,27);       % Number cloud drops deduced from optical thickness; number
aer1=xyz(:,28);     % Aerosol number, size>.1 micron and size<.2 micron; number
aer2=xyz(:,29);     % Aerosol number, size>.3 micron and size<.5 micron; number
aer3=xyz(:,30);     % Aerosol number, size>1 micron and size<5 micron; number
wvp=xyz(:,31);      % water vapor path; cm
qa=xyz(:,32);       % surface air specific humidity; g/kg
tlcl=xyz(:,33);     % lifting condensation level temperature; K
zlcl=xyz(:,34);     % lifting condensation level height; km
aert=xyz(:,35);     % total accumulation mode aerosols; number per cm^3
rh=xyz(:,36);       % surface relative humidity; %
pres=xyz(:,37);     % sea-level pressure; hPa
aerai=xyz(:,38);    % Aitken mode aerosol concentration; number per cm^3

% height-correct Z to align cloud tops
% height coordinate is middle of range gate
% triple height resolution to resolve transitions (e.g. top height) better in interpolation
height3=zeros(size(Z.height).*[3 1]);
height3(2:3:end-1)=Z.height;
height3(3:3:end-3)=(2*Z.height(1:end-1)+  Z.height(2:end))/3;
height3(4:3:end-2)=(  Z.height(1:end-1)+2*Z.height(2:end))/3;
height3(1)=        (4*Z.height(1      )-  Z.height(2    ))/3;
height3(end)=      (4*Z.height(end)    -  Z.height(end-1))/3;

ifct=isfinite(cloud_top_mean);
fifct=find(ifct);
mean_ct=mean(cloud_top_mean(ifct));
Z_cldtop_adj=NaN+zeros(size(Z.mean,1),length(height3));
% warning('off','MATLAB:interp1:NaNinY');
for i=1:sum(ifct) % slow, but no alternative
    Z_cldtop_adj(fifct(i),:)=interp1(Z.height'-cloud_top_mean(fifct(i))+mean_ct, Z.mean(fifct(i),:), height3);
end

[cmr,imax]=max(Z.mean,[],2);
[cmrsort,icmr]=sort(cmr);

cm=b2rcolormap(16);

% not height-adjusted
subplot(3,2,1,'align')
imagesc(time_yday,Z.height/1e3,Z.mean',[-40 20])
shading flat
%caxis([-40 20])
set(gca,'color',cm(1,:))
axis([316 336 0.1 2]); axis xy
title('VOCALS W-band reflectivity')
ylabel('height (km)')

subplot(3,2,2,'align')
imagesc(1:length(icmr),Z.height/1e3,Z.mean(icmr,:)',[-40 20])
shading flat
%caxis([-40 20])
set(gca,'color',cm(1,:))
axis([0 length(icmr) 0.1 2]); axis xy
title('sort by column max reflectivity')
hold on
plot((Z.height(imax(icmr))+25*rand(size(icmr))-12.5)/1e3,'m.','markersize',1)

% height-adjusted
subplot(3,2,3,'align')
imagesc(time_yday,height3/1e3,Z_cldtop_adj',[-40 20])
shading flat
%caxis([-40 20])
set(gca,'color',cm(1,:))
axis([316 336 0.1 1.6]); axis xy
ylabel('adjust all to mean cloud top height')
xlabel('yearday')

subplot(3,2,4,'align')
imagesc(1:length(icmr),height3/1e3,Z_cldtop_adj(icmr,:)',[-40 20])
shading flat
%caxis([-40 20])
set(gca,'color',cm(1,:))
axis([0 length(icmr) 0.1 1.6]); axis xy
hold on
plot((Z.height(imax(icmr))+25*rand(size(icmr))-12.5-cloud_top_mean(icmr)+mean_ct)/1e3,'m.','markersize',1)

hc=colorbar('southoutside');
set(hc,'position',get(hc,'position')+[0 -.07 0 0])
set(get(hc,'xlabel'),'string','dBZ')

orient landscape
set(gcf,'color','w','inverthardcopy','off')
print('-dpng',[way_proc_images_wband 'preproc_svd_wband_1min.png'])
% delete hires images
print('-depsc',[way_proc_images_wband 'preproc_svd_wband_1min.eps'])


