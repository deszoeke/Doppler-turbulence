% proc_wband_cloudtop.m
% 2009-11-09 :: VOCALS 2008 :: Simon de Szoeke
%
% Process the W-band radar Doppler velocity data for turbulence spectra.
%
% 2009-09-11 Update to include height correction for antenna zenith angle, when available. h=range*cos(zen_angle)

%cd /Users/sdeszoek/Data/cruises/VOCALS_2008/RHB/Scientific_analysis/programs/VOCALS2008_programs/wband
%addpath('~/Data/cruises/VOCALS_2008/RHB/Scientific_analysis/programs/VOCALS2008_programs/');
run('~/Data/cruises/VOCALS_2008/RHB/Scientific_analysis/programs/VOCALS2008_programs/read_parameters');

dodissipation=false;

% All radar data files
momentfile=dir([way_raw_data_wband '*MMCRMom.nc']);
year=2008;
base_time=NaN+zeros(length(momentfile),1);
base_time_offset=base_time;
ntimes=base_time;
time_interval=10; % minutes

% Memory strategy: 
% Cycle reflectivity data through 1 ~2100 member array.
% Write 1 file of 1-min cloud top for each MMCR file.
% Later concatenate files and average to 10-minute cloud top.

% get sizes to allocate Z array
for fi=1:length(momentfile)
    filename=[way_raw_data_wband momentfile(fi).name];
    A=nc_getvarinfo(filename,'time_offset');
    ntimes(fi)=A.Size(1);
end
%height=nc_varget(filename,'Heights',[0 0],[1 -1]); % actually range
height=ncread(filename,'Heights',[1 1],[Inf 1]); % actually range
h=height(height<2e3);                              %
range_gate=h(2)-h(1); % 24.9576 m
Z=NaN+zeros(max(ntimes),sum(height<2e3)); % preallocate, to be recycled

% After best recalibration, Ken Moran says to subtract 1.72 dB to all recorded values.
% reviever gain is 1.72 higher than previously estimated, which reduces dBZ
% and dBmW signals. Note recalibration noise source has an error of +-0.6 dB.
dB_offset=-1.72; % dB
% define a noise level for the reflectivity based on Ken Moran's manual
radar_const=19.613; % dB_
min_detectable_signal=-113.304+dB_offset; % dBmW
analog_noise_signal=-120.74; % dBmW peak diagnosed by Simon
digital_noise_signal=-115.4;  % dBmW
noise_margin=3.5; % dB
adhocthreshold=-43+dB_offset;
dBZnoiseavg=20*log10(h)+radar_const+analog_noise_signal;
dBZnoisefloor=20*log10(h)+radar_const+digital_noise_signal+noise_margin;
noisefloor=max(adhocthreshold,dBZnoisefloor);

% turbulence parameters
C=0.65; % nondimensional constant for Kolmogorov spectrum E_zz(k_x)
% sampling parameters
fsampl=3.5; % Hz
fnyq=fsampl/2;
band=fnyq./[20 3]; % frequency band over which to compute epsilon

% setup spectrum objects
% window design:
% N is 2100 for a 10-minute interval
nwindow=420; % 2 minute
noverlap=0.5*nwindow;
nws=9;       % 9 complete windows in the 10-minute interval
nfft=512;
% use spectrum objects to compute the spectrum
hwelch=spectrum.welch('Hann',nwindow,100*noverlap/nwindow); % --deprecated

% loop files
%starter=find(base_time_yday>=310,1,'first');
%starter=460;
starter=583; % first time any nonmissing data results in calculations
for fi=starter:length(momentfile) % fi>=2;
    %%prflname=[way_raw_data_wband momentfile(fi-1).name]; % previous file
    filename=[way_raw_data_wband momentfile(fi).name];
    fprintf(1,'%s ',momentfile(fi).name)
    
    % found no need to read the end of the previous file to get an integral minute
    % base_time in seconds since 1970-1-1 00:00:00
    %base_time(fi)=nc_varget(filename,'base_time');
    base_time(fi)=ncread(filename,'base_time');

    base_time_offset(fi)=base_time(fi)-base_time(starter);
    % base_time_mld in matlab datenumber
    base_time_mld=double(base_time(fi))/86400 + datenum(1970,1,1,0,0,0);
    base_time_yday=base_time_mld-datenum(year,1,0,0,0,0);
    % time_offset in seconds
    %time_offset=nc_varget(filename,'time_offset');
    time_offset=ncread(filename,'time_offset');
    % seconds since base_time(starter), for concatentation
    time_offset_cat=base_time_offset(fi)+time_offset;
    time_offset_bin=floor(base_time(fi)/60)*60-base_time(fi) + (0:time_interval*60:ceil(time_offset(end)/60)*60)';  

    % reflectivity matrix
    %Z(1:ntimes(fi),:)=nc_varget_lapxm(filename,'Reflectivity',[0 0],[-1 length(h)])+dB_offset;
    Z(1:ntimes(fi),:)=ncread(filename,'Reflectivity',[1 1],[length(h) Inf])'+dB_offset;
    %+dB_offset 2013-04-28 SPdeS (1.72 dB more liberal cloud detection was used for the 2011 paper);
    %vel(1:ntimes(fi),:)=nc_varget_lapxm(filename,'MeanDopplerVelocity',[0 0],[-1 length(h)]);
    vel(1:ntimes(fi),:)=ncread(filename,'MeanDopplerVelocity',[1 1],[length(h) Inf])';
    
    % Use a reflectivity criterion noisefloor defined above that selects clouds
    % and eliminates (almost) all noise returns.
    iscloud=Z(1:ntimes(fi),:)>repmat((noisefloor>0 & h<2e3 & h>400)',[ntimes(fi) 1]); % Z criterion for cloud
    
    % Read Kongsberg motion compensation file
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
        ifk=zeros(length(time_offset),1);
    elseif length(kongfile)>1
        kongflag=2;
        fkongerr=fopen([way_raw_data_wband 'motion_adjT/kongsberg_log.txt'],'a+');
        fprintf(fkongerr,['Multiple files found: ' kongfile(:).name]);
        fclose(fkongerr);
        wdrop=-vel;
        ifk=zeros(length(time_offset),1);
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
        theta=asin(sintheta); % zenith angle (radians)
        costheta=sqrt(1-quad);
        %height=costheta * range;

        wdrop=-vel(1:length(time_offset),:);
        % ship heave-compensated drop vertical velocity, when heave is
        % available from the Kongsberg
        ifk=isfinite(kongw4radar);
        wdrop(ifk,:)=-vel(ifk,:)-repmat(kongw4radar(ifk),[1,length(h)]);
    end

    % There are no missing values in w, but w is all noise for clear skies.
%     % interpolate w in clear skies
%     w(~iscloud)=interp1(find(iscloud),w(iscloud),find(~iscloud));
     
    % setup 10 minute bins
    num_returns=NaN+zeros(length(time_offset_bin)-1,1);
    num_heave_crct=NaN+zeros(length(time_offset_bin)-1,1);
    num_cloud=NaN+zeros(length(time_offset_bin)-1,length(h));
    [epsilon,std_eps,wbar,ww,www]=deal(NaN+zeros(length(time_offset_bin)-1,length(h)));
    %loop 10-min bins
    for i=1:length(time_offset_bin)-1
        in_bin=time_offset>=time_offset_bin(i) & time_offset<time_offset_bin(i+1);
        num_returns(i)=sum(in_bin); % number of columns in time interval
        num_cloud(i,:)=sum(iscloud(in_bin,:));
        num_heave_crct(i)=sum(ifk(in_bin));
        iz=num_cloud(i,:)>0.8*num_returns(i) & num_returns(i)>0.8*time_interval*60*fsampl; % vertical logical true where there are clouds
        U=repmat(10,[1 length(h)]); % 10 m/s for now, replace with HRDL velocity later
        %U=nanmean(u_hrdl(in_bin_hrdl,z_wband);
        if sum(iz)>0
            for iiz=find(iz)
                wfill=wdrop(in_bin,iiz); % nothing done to fill noncloud returns, noise remains!
                % compute 10-min mean, variance and triple product of w velocity
                wt=wfill(iscloud(in_bin,iiz));
                wbar(i,iiz)=mean(wt);
                wp=wt-wbar(i,iiz);
                ww(i,iiz)=mean(wp.*wp);
                www(i,iiz)=mean(wp.*wp.*wp);
                % compute dissipation
                if dodissipation
                    hpsd=psd(hwelch,wfill,'NFFT',nfft,'Fs',fsampl); % hwelch, parameters set above --deprecated
                    inband=hpsd.frequencies>band(1) & hpsd.frequencies<=band(2);
                    f=hpsd.frequencies(inband);
                    Ef=hpsd.data(inband,:); % spectral power per frequency
                    % k=2*pi*hpsd.frequencies/U;
                    epsilon(i,iiz)=2*pi*C^(-3/2)*mean( (f.^(5/2)/U(iiz)).*Ef.^(3/2) );
                    std_eps(i,iiz)=2*pi*C^(-3/2)*std( (f.^(5/2)/U(iiz)).*Ef.^(3/2) );
                end
            end
        end
    end % loop 10-min bins
    
    time_aggreg=base_time_offset(fi)+time_offset_bin(1:end-1); % seconds
   
    % write 10-minute turbulence dissipation file
    if dodissipation
        tfile=[way_proc_data_wband 'turbulence/' momentfile(fi).name(1:11) 'ResolvedTurb.txt'];
        ff=fopen(tfile,'w');
        for i=1:length(time_aggreg)
            fprintf(ff,'%7.0f\t',time_aggreg(i));
            fprintf(ff,'%6i\t',num_returns(i));
            fprintf(ff,'%6i\t',num_heave_crct(i));
            fprintf(ff,'%6i\t',num_cloud(i,:));
            fprintf(ff,'%6.4e\t',epsilon(i,:));
            fprintf(ff,'%6.4e\t',std_eps(i,:));
            fprintf(ff,'%6i\t',wbar(i,:));
            fprintf(ff,'%6.4e\t',ww(i,:));
            fprintf(ff,'%6.4e\t',www(i,:));
            fprintf(ff,'\n');
        end
        fclose(ff);
    end
    
    %
    tfile=[way_proc_data_wband 'turbulence/' momentfile(fi).name(1:11) 'w_ww_www_bar.txt'];
    ff=fopen(tfile,'w');
    for i=1:length(time_aggreg)
        fprintf(ff,'%6i\t',wbar(i,:));
        fprintf(ff,'%6.4e\t',ww(i,:));
        fprintf(ff,'%6.4e\t',www(i,:));
        fprintf(ff,'\n');
    end
    fclose(ff);

    
end % loop files

if ~exist([way_proc_data_wband 'turbulence/ResolvedTurb_10min_2008310-336.txt'],'file')
    filecat([way_proc_data_wband 'turbulence/2008*ResolvedTurb.txt'],...
            [way_proc_data_wband 'turbulence/ResolvedTurb_10min_2008310-336.txt'])
end

if ~exist([way_proc_data_wband 'turbulence/w_ww_www_10min_2008310-336.txt'],'file')
    filecat([way_proc_data_wband 'turbulence/2008*w_ww_www_bar.txt'],...
            [way_proc_data_wband 'turbulence/w_ww_www_10min_2008310-336.txt'])
end

%load the whole record of turbulence dissipation data
nh=length(h);
A=load([way_proc_data_wband 'turbulence/ResolvedTurb_10min_2008310-336.txt']);
% missing value of -999 may be supplied (currently is not)
i=1;
time=A(:,i); i=i+1;
num_returns=A(:,i); i=i+1;
num_heave=A(:,i); i=i+1;
num_cloud=A(:,i:i+nh-1); i=i+nh;
dissipation=A(:,i:i+nh-1); i=i+nh;
std_dissipation=A(:,i:i+nh-1); i=i+nh;
