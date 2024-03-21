% proc_wband_1min_stat_v6.m
% VOCALS 2008 :: 2009-07-07 :: Simon de Szoeke
%
% Compute 1-minute statistics of reflectivity and motion-compensated Doppler
% velocity from W-band moment files. 
%
% Statistics are vertically resolved. Height ranges (125-2000 m) are chosen
% for VOCALS 2008 clouds.
%
% Compute mean, min, max, variance, and skewness
% each minute for reflectivity and vertical velocity (+up)
%
% Histograms of reflectivity and w are computed as a function of height
% each minute.
%
% Output serial unformatted binary files, a fast intermediate format
% to be concatenated and compiled in a netcdf file by compile_wband_1min_stat.m.
%
% Version 2: Update calculation of the moments, comment out the histograms.
% Add housekeeping for number of observations.
%
% V.3 Subtract noise from received relectivity to get meteorological
% reflectivity. Record radar antenna zenith angle.
% V.4 Add Doppler spectral width statistics.
% V.5 Adjust by -1.72 dB from September calibration by Ken Moran.
% V.6 Average reflectivity Ze in linear mm^6/m^3 space 2010-07-13 (better for microphysical retrievals)
% 6/7 condition width on reflectivity>noise; used _onlywidth to get intermediate data files

% Dimensions:
% time( )                               serial time (s), record dimension
% range(75)     134:1982                center of range gates (m)
% Zbins(38)     [-Inf -40:2:30 Inf]     reflectivity histogram bin edges (dBZ)
% wbins(53)     [-Inf -5:0.2:5 Inf]     velocity histogram bins edges (m/s)
% Dwbins(73)    [-Inf  0:0.1:7 Inf];    Doppler vel width hist bin edges (m/s?)
read_parameters;

% dBZ and w bin edges
Zbins=[-inf -40:2:30 inf];
wbins=[-inf -5:0.2:5 inf];
Dwbins=[-inf 0:0.1:7 inf];

% All radar data files
momentfile=dir([way_raw_data_wband '*MMCRMom.nc']);
year=2008;
base_time=NaN+zeros(length(momentfile),1);
base_time_offset=base_time;
ntimes=base_time;

% get sizes to allocate refl array
for fi=1:length(momentfile)
    filename=[way_raw_data_wband momentfile(fi).name];
    A=nc_getvarinfo(filename,'time_offset');
    ntimes(fi)=A.Size(1);
end
range=nc_varget_lapxm(filename,'Heights',[0 0],[1 -1]);
hi=range>125 & range<2e3;
hi0=find(hi,1,'first');
h=range(hi);
refl=NaN+zeros(max(ntimes),length(h)); % preallocate, to be recycled
vel=refl;
wid=refl;

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
adhocthreshold=-43;
dBZnoiseavg=20*log10(h)+radar_const+analog_noise_signal;
dBZnoisefloor=20*log10(h)+radar_const+digital_noise_signal+noise_margin;
noisefloor=max(adhocthreshold,dBZnoisefloor);

% I find that the minimum detectable signal power of
% -113.304+dB_offset dBmW is pretty
% close to the 2nd (digital) noise peak, ~-111.6+dB_offset for day 326 hour 4.
% Do not subtract noisefloor from the received reflectivity to get the
% meteorological reflectivity, subtract the average noise. Take
% analog_noise_signal=-117.3 dBmW as the average noise power.
% This does not account for the finite width of the noise, which broadens the
% reflectivity signals.

% loop files
starter=460;
for fi=starter:length(momentfile) % fi>=2;
    filename=[way_raw_data_wband momentfile(fi).name];
    
    %%npre=300; % read a bit of the previous file to get an integral minute
    % base_time in seconds since 1970-1-1 00:00:00
    base_time(fi)=nc_varget(filename,'base_time');
    base_time_offset(fi)=base_time(fi)-base_time(starter);
    % base_time_mld in matlab datenumber and yearday
    base_time_mld=double(base_time(fi))/86400 + datenum(1970,1,1,0,0,0);
    base_time_yday=base_time_mld-datenum(year,1,0,0,0,0);
    % time_offset in seconds
    time_offset=nc_varget(filename,'time_offset');
    % seconds since base_time(starter), for concatentation
    time_offset_cat=base_time_offset(fi)+time_offset;

    % Read radar reflectivity and Doppler w
    refl(1:ntimes(fi),:)=nc_varget_lapxm(filename,'Reflectivity',[0 hi0-1],[-1 length(h)])+dB_offset;  % Reflectivity, dBZ
    vel(1:ntimes(fi),:)=nc_varget_lapxm(filename,'MeanDopplerVelocity',[0 hi0-1],[-1 length(h)]);      % Doppler velocity, m/s
    wid(1:ntimes(fi),:)=nc_varget_lapxm(filename,'SpectralWidth',[0 hi0-1],[-1 length(h)]);            % Doppler spectral width, m/s

    Zrec=10.^(refl./10); % retreived Ze, mm^6/m^3
    % subtract off average noise
    % This biases result! Instead subtract off average noise from averages below.
%     Zmet=Zrec - repmat(10.^(dBZnoiseavg./10),[size(refl,1) 1]);
%     reflmet=real(10*log10(Zmet)); % dBZ from meteorological backscatter
%     reflmet(Zmet<0)=NaN; % biases results high!
    
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
        fprintf(fkongerr,['No Kongsberg file in hour:' ddd ' ' hh '\n']);
        fclose(fkongerr);
        wdrop=-vel;
    elseif length(kongfile)>1
        kongflag=2;
        fkongerr=fopen([way_raw_data_wband 'motion_adjT/kongsberg_log.txt'],'a+');
        fprintf(fkongerr,['Multiple files found: ' kongfile(:).name '\n']);
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
        theta=asin(sintheta); % zenith angle (radians)
        costheta=sqrt(1-quad);
        %height=costheta * range;

        % motion-compensated drop vertical velocity +up
        wdrop=-vel(1:length(time_offset),:)-repmat(kongw4radar,[1,length(h)]);
    end
    
    % 1-minute bins
    time_offset_bin=floor(base_time(fi)/60)*60-base_time(fi) + (0:60:ceil(time_offset(end)/60)*60)';
    Z.cfad=NaN+zeros(length(Zbins),length(h),length(time_offset_bin)-1);
    w.cfad=NaN+zeros(length(wbins),length(h),length(time_offset_bin)-1);
    Dw.cfad=NaN+zeros(length(Dwbins),length(h),length(time_offset_bin)-1);
    Z.mean=NaN+zeros(length(time_offset_bin)-1,length(h));
    Z.max=Z.mean;
    Z.min=Z.mean;
    Z.var=Z.mean;
    Z.skew=Z.mean;
    w.mean=NaN+zeros(length(time_offset_bin)-1,length(h));
    w.max=w.mean;
    w.min=w.mean;
    w.var=w.mean;
    w.skew=w.mean;
    Dw.mean=NaN+zeros(length(time_offset_bin)-1,length(h));
    Dw.max=Dw.mean;
    Dw.min=Dw.mean;
    Dw.var=Dw.mean;
    Dw.skew=Dw.mean;
    zen_angle.mean=NaN+zeros(length(time_offset_bin)-1,1);
    zen_angle.std =zen_angle.mean;
    zen_angle.cos =zen_angle.mean;
    zen_angle.sin =zen_angle.mean;
    %loop 1-minute bins
    for i=1:length(time_offset_bin)-1
        in_bin=time_offset>=time_offset_bin(i) & time_offset<time_offset_bin(i+1);       
        if sum(in_bin)>0
            % mask out velocities that have sub-noise reflectivity
            mask=refl(in_bin,:)>repmat(noisefloor,[sum(in_bin) 1]);
            wdropmask=wdrop(in_bin,:);
            wdropmask(~mask)=NaN;
            % condition width on reflectivity>noise too 7/5/2011
            widmask=wid(in_bin,:);
            widmask(~mask)=NaN;
            
            % counts in bins for CFADs
            Z.cfad(:,:,i)=histc(refl(in_bin,:),Zbins); % dBZ
            w.cfad(:,:,i)=histc(wdropmask,wbins);
            Dw.cfad(:,:,i)=histc(widmask,Dwbins);
            
            % statistics over time: mean max min variance skewness
            % V.6 2010-07-13 average Z in linear space rather than dBZ, then convert to dBZ
            Z.mean(i,:)=10*log10( nanmean(Zrec(in_bin,:)) );
            Z.min(i,:) =10*log10( min(Zrec(in_bin,:)) );
            Z.max(i,:) =10*log10( max(Zrec(in_bin,:)) );
            Z.var(i,:) =10*log10( nanvar(Zrec(in_bin,:),1) ); % second-moment
            Z.skew(i,:)=10*log10( skewness(Zrec(in_bin,:)) );
            
            w.mean(i,:)=nanmean(wdropmask);
            w.min(i,:) =min(wdropmask);
            w.max(i,:) =max(wdropmask);
            w.var(i,:) =nanvar(wdropmask,1); % second-moment
            w.skew(i,:)=skewness(wdropmask);
            
            Dw.mean(i,:)=nanmean(wid(in_bin,:));
            Dw.min(i,:) =min(wid(in_bin,:));
            Dw.max(i,:) =max(wid(in_bin,:));
            Dw.var(i,:) =nanvar(wid(in_bin,:),1); % second-moment
            Dw.skew(i,:)=skewness(wid(in_bin,:));
            
            if ~kongflag
                zen_angle.mean(i)=nanmean(theta(in_bin));
                zen_angle.std (i)=nanstd (theta(in_bin));
                zen_angle.cos (i)=nanmean(costheta(in_bin));
                zen_angle.sin (i)=nanmean(sintheta(in_bin));
            else
                zen_angle.mean(i)=NaN;
                zen_angle.std (i)=NaN;
                zen_angle.cos (i)=NaN;
                zen_angle.sin (i)=NaN;
            end
        end
    end
    
    % write 1-minute data file for each W-band Mom file
    time_aggreg=base_time_offset(fi)+time_offset_bin(1:end-1); % seconds
    %sprintf('%5.0f\t%6.4e\t%9.4f\t%3d\t%1d\n',[time_aggreg, cloud_fraction, mean_top_bin, double(num_returns), double(flag)]')
    
    writetimeserially=true;
    
    % Save 1min_stat for each Mom file
    statfile=[way_proc_data_wband '1min_stat/' momentfile(fi).name(1:11) '_1min_stat.dat'];
    fs=fopen(statfile,'w','ieee-le');
    if ~writetimeserially    % write variables serially
        fwrite(fs,time_aggreg,'int32');  % 1 integer 32-bit
        fwrite(fs,[zen_angle.mean zen_angle.std zen_angle.cos zen_angle.sin],'float32'); % 4 floats 32-bit
        fwrite(fs,[Z.mean Z.min Z.max Z.var Z.skew w.mean w.min w.max w.var w.skew Dw.mean Dw.min Dw.max Dw.var Dw.skew],'float64');  % 10 floats 64-bit
    else % write time serially - suitable for concatenating files along time dimension
        A=[Z.mean Z.min Z.max Z.var Z.skew w.mean w.min w.max w.var w.skew Dw.mean Dw.min Dw.max Dw.var Dw.skew];
        for i=1:length(time_offset_bin)-1
            fwrite(fs,time_aggreg(i),'int32');  % 1 integer 32-bit
%             fwrite(fs,kongflag,'int8'); % 0 good, 1 no file, 2 multiple files
            fwrite(fs,[zen_angle.mean(i) zen_angle.std(i) zen_angle.cos(i) zen_angle.sin(i)],'float32'); % 4 floats 32-bit
            fwrite(fs,A(i,:),'float64');        % 10 floats 64-bit
        end
    end
    fclose(fs);
    
    % Save CFAD for each Mom file
    cfadfile=[way_proc_data_wband '1min_stat/' momentfile(fi).name(1:11) '_cfad.dat'];
    fc=fopen(cfadfile,'w','ieee-le');
    if ~writetimeserially % write variables serially
        fwrite(fc,Z.cfad,'uint8');
        fwrite(fc,w.cfad,'uint8');
        fwrite(fc,Dw.cfad,'uint8');
    else % write time serially - suitable for concatenating files along time dimension
        for i=1:length(time_offset_bin)-1
            fwrite(fc,Z.cfad(:,:,i),'uint8');
            fwrite(fc,w.cfad(:,:,i),'uint8');
            fwrite(fc,Dw.cfad(:,:,i),'uint8');
        end
    end
    fclose(fc);

end % loop Mom files (fi)

%save range variable (formerly called height, 2009 September 10)
rangefile=[way_proc_data_wband '1min_stat/range.txt'];
fh=fopen(rangefile,'w');
fprintf(fh,'%6.3f\n',h);
fclose(fh);
