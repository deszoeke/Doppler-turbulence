% Retrieve microphyscial and macrophysical quantities for drizzle and clouds
% in drizzling stratocumulus of the southeastern tropical Pacific for 
% VOCALS 2008.
% Inputs:  Microwave radiometer LWP, W-band radar cloud top height,
%          ceilometer cloud base height, W-band radar Doppler moments
%          (Z, velocity, velocity spectral width), Kongsberg radar
%          motion, HRDL lidar horizontal velocity.
% Outputs: Detect cloud or drizzle from radar moments.
%          Drop number concentration, liquid water content,
%          drop fall velocity, lognormal modal radius
%
% VOCALS 2008 :: 2009-09-23 :: Simon de Szoeke

read_parameters;

% Define dimensions
% vertical coordinate (center of ranges, m)
% h=load([way_proc_data_wband '1min_stat/height.txt']);
% time_aggreg in seconds since base_time
base_time=1225916508;
% time_yday=datenum(0,0,0,0,0,base_time+time_aggreg)+datenum(1970,1,1,0,0,0)-datenum(2008,1,0,0,0,0); %yearday

load([way_proc_data_wband '1min_stat/Z_1min.mat']);

nheight=length(Z.height); % 75
nmin=60;
nhour=24;

% load 1-minute cloud fraction and top height data derived from W-band (Simon de Szoeke, Oregon State U.)
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
time_yday=time_top_yday;

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

% synchronize time bases
% 1-minute time base
temp=time_yday(1)+0:1/(60*24):time_yday(end);  % even 1-min time_yday without gaps
time1yday=NaN+zeros(10*ceil(length(temp)/10),1); % pad w/ nans to get length multiple of 10
time1yday(1:length(temp))=temp;
n1m=length(time1yday);
% 10-minute time base
time10yday=time1yday(1:10:end); % even 10-min
n10m=length(time10yday);

% synchronize 1-minute data
% time1yday(loc(ii)) approx. equals time_yday(ii)

% synchronize 1-min wband radar cloud top to cloudtop1(time1yday)
[ii,loc]=ismember(round(time_yday*60*24),round(time1yday*60*24));
cloudtop1=NaN+time1yday;
cloudtop1(loc(ii))=cloud_top_valmean(ii);

% % synchronize 10-min W-band cloud top height cth10 to cloudtop10(time10yday)
% [ii,loc]=ismember(round(cth10.yday*60*24/10),round(time10yday*60*24/10));
% cloudtop10=NaN+time10yday;
% cloudtop10(loc(ii))=cth10.height(ii);

%synchronize ceilo cloud base height clb10 to cloudbase10(time10yday)
[ii,loc]=ismember(round(clb10.yday*60*24/10),round(time10yday*60*24/10));
cloudbase10=NaN+time10yday;
cloudbase10(loc(ii))=clb10.height(ii);

% cloudbase10(cloudbase10>cloudtop10)=NaN;
% cloudtop10(cloudbase10>cloudtop10)=NaN;

cloudbase1=ave_large2small(time1yday,time10yday,cloudbase10);
cloudtop1(cloudtop1<cloudbase1)=NaN;
cloudbase1(cloudtop1<cloudbase1)=NaN;

% synchronize LWP and Ztopmax10(time10yday)
[ii,loc]=ismember(round(lwp.yday*60*24/10),round(time10yday*60*24/10));
% time10yday(loc(ii)) approx equals lwp.yday(ii)
lwp10=NaN+cloudbase10;
adlwp10=lwp10;
lwp10(loc(ii))=lwp.lwp(ii);
adlwp10(loc(ii))=lwp.adlwp(ii);
lwp1=ave_large2small(time1yday,time10yday,cloudbase1);

% compute 1-minute volume mixing ratio at top from LWP and cloud thickness
Qctop=2*lwp1./(cloudtop1-cloudbase1);
% mass mixing ratio
rhosfc=1.19;
rhotop=rhosfc*exp(-cloudtop1/8); % kg m^-3, approx from 8 km scale height
qctop=Qctop./rhotop;

% After best recalibration, Ken Moran says to subtract 1.72 dB to all recorded values.
% reviever gain is 1.72 higher than previously estimated, which reduces actual dBZ
% and dBmW signals. Note recalibration noise source has an error of +-0.6 dB.
dB_offset=-1.72; % dB
% define a noise level for the reflectivity based on Ken Moran's manual
radar_const=19.613; % dB_
min_detectable_signal=-113.304+dB_offset; % dBmW
adhocthreshold=-43+dB_offset;
% I find that the minimum detectable signal power of -113.304 dBmW is pretty
% close to where I see the 2nd noise peak, ~-113.75 for day 326 hour 4.
%
% Take 113.304 dBmW as the noise power and subtract it range-corrected (dBZnoise)
% from the received reflectivity in linear units to get the meteorological reflectivity.
% This does not account for the finite width of the noise, which broadens the
% reflectivity signals. Also, still threshold the data below noisefloor.

% load 1 hour of Wband radar data
starter=460;
year=2008;
yday=325;
hour=12;
fi=1;
    yyyy=sprintf('%04d',year);
    ddd= sprintf('%03i',floor(yday));
    hh=  sprintf('%02d',floor(hour));

    file=dir([way_raw_data_wband yyyy ddd hh '*MMCRMom.nc']);
    filename=[way_raw_data_wband file.name];
    % base_time in seconds since 1970-1-1 00:00:00
    base_time(fi)=nc_varget(filename,'base_time');
    %base_time_offset(fi)=base_time(fi)-base_time(starter);
    % base_time_mld in matlab datenumber and yearday
    base_time_mld=double(base_time(fi))/86400 + datenum(1970,1,1,0,0,0);
    base_time_yday=base_time_mld-datenum(year,1,0,0,0,0);
    % time_offset in seconds
    time_offset=nc_varget(filename,'time_offset');
    % seconds since base_time(starter), for concatentation
    %time_offset_cat=base_time_offset(fi)+time_offset;

    % height coordinate
    range=nc_varget_lapxm(filename,'Heights',[0 0],[1 -1]);
    hi=range>125 & range<2e3;
    hi0=find(hi,1,'first');
    h=range(hi);
    
    % range-dependent noise floor
    dBZnoise=20*log10(h)+radar_const+min_detectable_signal;
    noisefloor=max(adhocthreshold,dBZnoise);

    % Read radar reflectivity and Doppler w
    refl=nc_varget_lapxm(filename,'Reflectivity',[0 hi0-1],[-1 length(h)]) + ... % Reflectivity, dBZ
         +dB_offset;
    vel=nc_varget_lapxm(filename,'MeanDopplerVelocity',[0 hi0-1],[-1 length(h)]); % Doppler velocity, m/s
    wid=nc_varget_lapxm(filename,'SpectralWidth',[0 hi0-1],[-1 length(h)]);       % Doppler spectral width, m/s
    
    % subtract off mean noise
    Zrec=10.^(refl./10); % received Z [mm6/m3]
    Zmet=Zrec - repmat(10.^(dBZnoise./10),[size(refl,1) 1]);
    reflmet=real(10*log10(Zmet)); % dBZ from meteorological backscatter
    reflmet(Zmet<0)=NaN;
    
    % Read 1-hour Kongsberg motion compensation file
    yyyy=sprintf('%04d',year);
    ddd= sprintf('%03i',floor(base_time_yday));
%     hh=  sprintf('%02d',floor(mod(base_time_yday,1)*24)); % subject to rounding error!
    DV=datevec(base_time_mld);
    hh=  sprintf('%02d',DV(4));
    kongfile=dir([way_raw_data_wband 'motion_adjT/' yyyy ddd hh '*Kongsberg_adjT.txt']);
    if isempty(kongfile)
        kongflag=1;
        fkongerr=fopen([way_raw_data_wband 'motion_adjT/kongsberg_log.txt'],'a+');
        fprintf(fkongerr,['No Kongsberg file in hour:' ddd ' ' hh]);
        fclose(fkongerr);
        wdrop=-vel;
    elseif length(kongfile)>1
        kongflag=2;
        fkongerr=fopen([way_raw_data_wband 'motion_adjT/kongsberg_log.txt'],'a+');
        fprintf(fkongerr,['Multiple files found: ' kongfile(:).name]);
        fclose(fkongerr);
        wdrop=-vel;
    else
        kongflag=0;
        kongfilename=[way_raw_data_wband 'motion_adjT/' kongfile(1).name];
        [kongtime,kongw,kongpitch,kongroll]=read_kongsberg(kongfilename);
        
        % synchronize Kongsberg and radar time (yearday)
        [junk,ii,jj]=unique(kongtime);
        kongw4radar=interp1(kongtime(ii),kongw(ii)    ,base_time_yday+time_offset/86400);
        pitch      =interp1(kongtime(ii),kongpitch(ii),base_time_yday+time_offset/86400); %degrees
        roll       =interp1(kongtime(ii),kongroll(ii) ,base_time_yday+time_offset/86400); %      
        % pitch/roll housekeeping vertical coordinate
        quad=sin(pitch/180*pi).^2+sin(roll/180*pi).^2; % equals sin(theta)^2
        sintheta=sqrt(quad);
        costheta=sqrt(1-quad);
        theta=asin(sintheta); % zenith angle (radians)
        %height=costheta * range;

        % motion-compensated drop vertical velocity +up
        wdrop=-vel(1:length(time_offset),:)-repmat(kongw4radar,[1,length(h)]);
        
        % wind compensation for horiz. (u,v) when antenna tilts...?
    end
    % Doppler velocity (+toward) ~= fall speed (+down) ~= -wdrop (+up)

% Assume DRIZZLE with lognormal radius distribution and Doppler velocity=Fall speed
[sigx_driz, r0_driz, N_driz]=drizzle_lognorm_param(-wdrop,wid,Zmet*1e-18); % Z[mm6/mm3]-->[m3]
rhow=1e3; % liquid water density kg/m^3;
lwc_driz=pi/6*rhow*sqrt(N_driz.*Zmet*1e-18).*exp(-9/2*sigx_driz.^2)*1e3; % g/m^3
a=1.2e-4; % s
b=1.0e-5; % m
xi=reflmet-180+10*log10(pi*1e3*rhow/48)+10*6*log10(exp(1))*sigx_driz.^2; % ~(reflmet-132+26*sigx_driz.^2)
ql_16=(-wdrop*a+b).^(-3).*10.^(xi/10); % liquid water mixing ratio (Frisch et al. 1995)
Fq=-ql_16.*((-wdrop+b/a).*exp(-3*xi.^2)-b/a);

figure; dock; clf
subplot(3,1,1)
imagesc((time_offset-time_offset(1))/60,h/1e3,sigx_driz')
axis xy
caxis([0 3])
title({['Drizzle lognormal retrieval, VOCALS 2008, day ' ddd ' hour ' hh], 'modal width \sigma_{ln(r)}'})
colorbar
subplot(3,1,2)
imagesc((time_offset-time_offset(1))/60,h/1e3,r0_driz'*1e3)
axis xy
caxis([0 0.2]);
colorbar
title('modal radius (mm)')
subplot(3,1,3)
imagesc((time_offset-time_offset(1))/60,h/1e3,log10(max(0,N_driz)'*1e-6))
axis xy
%caxis([0 1000])
colorbar
title('drizzle drop number (log_{10} cm^{-1})')
%print('-dpng',[way_proc_images_wband 'drizzle_retr1_0_' yyyy ddd hh '.png'])

% Note color scales are influenced by outliers. Histogram may be more
% informative.
% From the plot, this retrieval often gives nonsensical results.