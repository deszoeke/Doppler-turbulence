% diurnal_CMR.m
% VOCALS 2008 :: 2010-06-18 :: Simon de Szoeke
%
% Joint diurnal-column max reflectivity histogram for 1-min VOCALS leg 2 W-band radar reflectivity.

% Requires data:
% cloudheight/CloudHeight_1min_2008310-336.txt from proc_wband_cloudtop.m
% 1min_stat/Z_1min.mat from compil_1min_stat.m
% ~/Data/Stratus/cloud_10min2008.txt VOCALS synthesis data

% preamble
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

% joint histogram of time of day and CMR
lon1=interp1(jd,lon,time_yday);
todh=mod(time_yday+lon1/360,1)*24; % time of day hour
cmr_edge=cmrsort(1:floor(length(cmr)/31):length(cmr));
cmr_edge(end)=max(cmr)+eps;
for hi=1:24
    for ri=1:length(cmr_edge)-1
    ii=todh>=hi-1 & todh<hi & cmr>=cmr_edge(ri) & cmr<cmr_edge(ri+1);
    cnt(hi,ri)=sum(ii);
    end
end

imagesc(-24:23,1:31,cnt([1:end 1:end],:)');
axis xy
[ucmr,ui]=unique(cmrsort(isfinite(cmrsort)));
ydbz=[-35:5:-20 -10:10:10];
ytk=interp1(ucmr,ui,ydbz);
set(gca,'tickdir','out','fontsize',14,'xtick',-24.5:6:23.5,'xticklabel',0:6:18,'ytick',sort([1:3:31])-.5,'yticklabel',ydbz)
set(gca,'yticklabel',round(cmr_edge([1:3:31])))
xlabel('local hour')
ylabel('CMR dBZ')
colorbar('southoutside')
% stuff taken off...
% [c,h]=contourf((0.5:31)*1071,-48.5:24.5,cnt([end 1:end 1:end 1:end 1],[1:end]));
% set(h,'edgecolor','none')
% hold on
% ax(1)=gca;
% ax(2)=axes('position',get(ax(1),'position')-[0 .05 0 0],'visible','off');
% set(ax(1:2),'xlim',[-48.5 24.5])
% plot(todh(icmr),1:length(cmr),'b.','markersize',1)
% plot(todh(icmr)-24,1:length(cmr),'b.','markersize',1)
% plot(todh(icmr)-48,1:length(cmr),'b.','markersize',1)
% plot(cmr(icmr),1:length(cmr),'r-','linewidth',1.2); % cumulative distribution
% ax(2)=axes('position',get(ax(1),'position').*[1 1 72/73 eps]-[0 .1 0 0],'xlim',[-48 24],'xtick',-48:6:24,'xticklabel',0:6:18);
% set(get(ax(2),'xlabel'),'string','local hour')
% ax(3)=axes('position',get(ax(1),'position').*[1 1 eps 30/31]-[.1 0 0 0],'ylim',[0 length(ui)],'ytick',ytk,'yticklabel',ydbz);
% set(get(ax(3),'ylabel'),'string','dBZ');

orient portrait
print('-depsc',[way_proc_images_wband 'diurnal_CMR.eps'])


