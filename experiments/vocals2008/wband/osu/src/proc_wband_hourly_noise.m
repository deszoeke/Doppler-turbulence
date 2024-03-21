% proc_wband_hourly_noise.m
% VOCALS 2008 :: 2009-10-19 :: Simon de Szoeke
% 
% Compute hourly noise power histograms for W-band radar.

read_parameters;

binrucnoise=-122:.02:-100;

% All radar data files
momentfile=dir([way_raw_data_wband '*MMCRMom.nc']);
year=2008;
base_time=NaN+zeros(length(momentfile),1);
base_time_offset=base_time;
ntimes=base_time;

% get sizes to allocate refl array
for fi=1:length(momentfile)
    filename=[way_raw_data_wband momentfile(fi).name];
end
range=nc_varget_lapxm(filename,'Heights',[0 0],[1 -1]);
hi=range>125 & range<2e3;
hi0=find(hi,1,'first');
h=range(hi);

% After best recalibration, Ken Moran says to subtract 1.72 dB to all recorded values.
% reviever gain is 1.72 higher than previously estimated, which reduces dBZ
% and dBmW signals. Note recalibration noise source has an error of +-0.6 dB.
dB_offset=-1.72; % dB
% define a noise level for the reflectivity based on Ken Moran's manual
radar_const=19.613; % dB_
min_detectable_signal=-113.304+dB_offset; % dBmW
analog_noise_signal=-120.74; % dBmW peak diagnosed by Simon
digital_noise_signal=-115.4;
noise_margin=3.5; % dB
adhocthreshold=-43+dB_offset;
dBZnoiseavg=20*log10(h)+radar_const+analog_noise_signal;
dBZnoisefloor=20*log10(h)+radar_const+min_detectable_signal+noise_margin;
noisefloor=max(adhocthreshold,dBZnoisefloor);

% file to fill with hourly histograms of near-noise power, dBmW
fidnoise=fopen([way_proc_data_wband 'noise_dBmW.dat'],'w','ieee-le');
fidnoisetime=fopen([way_proc_data_wband 'noise_time.dat'],'w','ieee-le');

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
    RCpower=nc_varget_lapxm(filename,'RangeCorrectedPower',[0 hi0-1],[-1 length(h)]) +dB_offset;            % Doppler spectral width, m/s
    
    % hourly noise & near-noise signal histograms
    % save CFAD of (range-UNcorrected range corrected) calibrated power, dBmW
    hh=repmat(h,[length(time_offset) 1]);
    rucpower=RCpower-20*log10(hh);
    countnoise=histc(rucpower(hh>100 & hh<=3000),binrucnoise); % store/save as uint16
    fwrite(fidnoise,countnoise,'uint16');
    fwrite(fidnoisetime,time_offset_cat(1),'uint32');

end % loop Mom files (fi)
fclose(fidnoise);
fclose(fidnoisetime);
