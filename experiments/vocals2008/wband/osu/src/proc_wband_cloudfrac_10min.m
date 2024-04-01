% proc_wband_cloudfrac_10min
% VOCALS 2008 :: 2009-10-21 :: Simon de Szoeke
%
% Average cloud fraction from W-band radar (1-minute resolution)
% to 10 minute resolution. Requires output from proc_wband_cloudfraction.m
% to be available first.

% preamble
read_parameters

% load the 1-min cloud fraction and data
A=load([way_proc_data_wband 'cloudfraction/CloudFraction_1min_2008310-336.txt']);
A(A==-999)=NaN; % NaN out missing value of -999
time=A(:,1);
cloud_frac1=A(:,2);
cloud_frac2=A(:,3);
cloud_frac3=A(:,4);
glitch_frac=A(:,5);

T=load([way_proc_data_wband 'cloudfraction/ClFrThresh_1min_2008310-336.txt']);
% same time base as A
cf=T(:,2:end);
thr=[0 0.5 1 2 4 8 12];

starter=1;
base_time(starter)=1225916508; % s since 1970-1-1 00:00
time_yday=datenum(0,0,0,0,0,base_time(starter)+time)+datenum(1970,1,1,0,0,0)-datenum(2008,1,0,0,0,0);

% AVERAGE to 10-minute (600 s) interval
sec1m=60;
sec10min=60*10;
% varaibles are named: name_units_resolution: sec=seconds, 1m=1 minute, 10m=10 minutes
time_sec_1m=base_time(starter)+time;
time0_10m_1m=floor((base_time(starter)+time)/sec10min);
time0_10m_10m=unique(time0_10m_1m);
time0_sec_10m=unique(time0_10m_1m)*sec10min; %starting time in seconds of 10m interval
time0_yday_10m=datenum(0,0,0,0,0,time0_sec_10m)+datenum(1970,1,1,0,0,0)-datenum(2008,1,0,0,0,0);
n10m=length(unique(time0_10m_10m));
cloud_frac1_10m=NaN+zeros(n10m,1);
cloud_frac2_10m=NaN+zeros(n10m,1);
cloud_frac3_10m=NaN+zeros(n10m,1);
glitch_frac_10m=NaN+zeros(n10m,1);
cf_10m=NaN+zeros(n10m,size(cf,2));
for ti=1:n10m
    iitime=time0_10m_1m==time0_10m_10m(ti);
    cloud_frac1_10m(ti)=nanmean(cloud_frac1(iitime)); % weights each minute the same
    cloud_frac2_10m(ti)=nanmean(cloud_frac2(iitime)); % number of cloudy returns
    cloud_frac3_10m(ti)=nanmean(cloud_frac3(iitime));
    glitch_frac_10m(ti)=nanmean(glitch_frac(iitime));
    cf_10m(ti,:)=nanmean(cf(iitime,:));
end

% plot to check
plot(time0_yday_10m,cloud_frac1_10m,'r.')
hold on
plot(time_yday,cloud_frac1,'-','markersize',3)

% save 10-min cloud fraction in a text file
fracfile=[way_proc_data_wband 'cloudfraction/VOCALS2008_WbandCloudFraction10min_1_0.txt'];
ff=fopen(fracfile,'w');
fprintf(ff,'%12.8f\t%10.5f\t%10.5f\t%10.5f\t%10.5f\n',...
    [time0_yday_10m, max(-999,cloud_frac1_10m), max(-999,cloud_frac2_10m),...
     max(-999,cloud_frac3_10m), max(-999,glitch_frac_10m)]');
fclose(ff);
% save other threshold cloud fractions in a text file
thrfile=[way_proc_data_wband 'cloudfraction/VOCALS2008_WbandClFrThresh10min.txt'];
ff=fopen(thrfile,'w');
fprintf(ff,'%12.8f\t%10.5f\t%10.5f\t%10.5f\t%10.5f\t%10.5f\t%10.5f\t%10.5f\n',...
    [time0_yday_10m, max(-999,cf_10m)]')
fclose(ff);
    
% to read 10-minute data
%fracfile=[way_proc_data_wband 'cloudfraction/CloudFraction_10min.txt'];
missing_value=-999;
ff=fopen(fracfile);
B=fscanf(ff,'%f',[7 3234])';
B(B==missing_value)=NaN;
yday=B(:,1);
cldfrac1=B(:,2);
cldfrac2=B(:,3);
cldfrac3=B(:,4);
glitchfr=B(:,5);
fclose(ff);

% save 10-minute cloud fraction in a netcdf file
ncpath=[way_proc_data_wband 'cloudfraction/']; % note netcdflib does not accept ~/ as a path
%ncf=[ncpath 'VOCALS2008_WbandCloudFraction10min_1_0.nc'];
ncf=[ncpath 'VOCALS2008_WbandCloudFraction10min_1_0a.nc']; % v1.0a
nc_create_empty( ncf );
nc_add_dimension( ncf, 'time', 0 ) % record dimension
nc_add_dimension( ncf, 'threshold', size(cf_10m,2)) % v1.0a

varstruct.Name = 'yday';
varstruct.Nctype = nc_float;
varstruct.Dimension = {'time'}; % { cell } needed
nc_addvar ( ncf , varstruct )
% more floats...
varstruct.Name = 'cloud_fraction1'; nc_addvar ( ncf , varstruct );
varstruct.Name = 'cloud_fraction2'; nc_addvar ( ncf , varstruct );
varstruct.Name = 'cloud_fraction3'; nc_addvar ( ncf , varstruct );
varstruct.Name = 'glitch_fraction'; nc_addvar ( ncf , varstruct );
varstruct.Name = 'cloud_fraction_threshold';                        % v1.0a
varstruct.Dimension = {'time','threshold'};                         % v1.0a
nc_addvar ( ncf , varstruct );                                      % v1.0a
varstruct.Name = 'reflectivity_threshold';                          % v1.0a
varstruct.Dimension = {'threshold'};                                % v1.0a
nc_addvar( ncf , varstruct );                                       % v1.0a

% variable attributes
nc_attput ( ncf , 'yday', 'long_name', 'yearday since 2007-December-31 00:00 (1.0 for 2008 January 1 00:00)' );
nc_attput ( ncf , 'yday', 'units', 'day' );
nc_attput ( ncf , 'yday', 'description', 'beginning of time interval');
nc_attput ( ncf , 'reflectivity_threshold', 'long_name', 'reflectivity threshold for cloud detection' ); % v1.0a
nc_attput ( ncf , 'reflectivity_threshold', 'units', 'dBZ' );                                            % v1.0a
nc_attput ( ncf , 'reflectivity_threshold', 'description', 'dBZ above noise margin');                    % v1.0a

nc_attput ( ncf , 'cloud_fraction1', 'long_name', 'W-band cloud fraction 1');
nc_attput ( ncf , 'cloud_fraction1', 'description', 'fraction of columns with >=1 return gate above noise reflectivity');
nc_attput ( ncf , 'cloud_fraction1', 'units', 'fraction');
nc_attput ( ncf , 'cloud_fraction1', 'valid_range', [0 1]);
%nc_attput ( ncf , 'cloud_fraction1', 'missing_value', missing_value)
nc_attput ( ncf , 'cloud_fraction2', 'long_name', 'W-band cloud fraction 2');
nc_attput ( ncf , 'cloud_fraction2', 'description', 'fraction of columns with >=2 vertically contiguous return gates above noise reflectivity');
nc_attput ( ncf , 'cloud_fraction2', 'units', 'fraction');
nc_attput ( ncf , 'cloud_fraction2', 'valid_range', [0 1]);
%nc_attput ( ncf , 'cloud_fraction2', 'missing_value', missing_value)
nc_attput ( ncf , 'cloud_fraction3', 'long_name', 'W-band cloud fraction 3');
nc_attput ( ncf , 'cloud_fraction3', 'description', 'fraction of columns with >=3 vertically contiguous return gates above noise reflectivity');
nc_attput ( ncf , 'cloud_fraction3', 'units', 'fraction');
nc_attput ( ncf , 'cloud_fraction3', 'valid_range', [0 1]);
%nc_attput ( ncf , 'cloud_fraction3', 'missing_value', missing_value)
nc_attput ( ncf , 'glitch_fraction', 'long_name', 'W-band cloud fraction 1');
nc_attput ( ncf , 'glitch_fraction', 'description', 'fraction of columns with exactly 1 or 2 return gates above noise reflectivity');
nc_attput ( ncf , 'glitch_fraction', 'units', 'fraction');
nc_attput ( ncf , 'glitch_fraction', 'valid_range', [0 1]);
%nc_attput ( ncf , 'glitch_fraction', 'missing_value', missing_value)
nc_attput ( ncf , 'cloud_fraction_threshold', 'long_name', 'W-band cloud fraction with different reflectivity thresholds'); % v1.0a
nc_attput ( ncf , 'cloud_fraction_threshold', 'description', 'no vertical contiguousness condition'); % v1.0a
nc_attput ( ncf , 'cloud_fraction_threshold', 'units', 'fraction'); % v1.0a
nc_attput ( ncf , 'cloud_fraction_threshold', 'valid_range', [0 1]); % v1.0a

% write the data
nc_varput ( ncf, 'yday'                  , time0_yday_10m     );
nc_varput ( ncf, 'reflectivity_threshold', thr'    );         % v1.0a

nc_varput ( ncf, 'cloud_fraction1'   , cloud_frac1_10m    );
nc_varput ( ncf, 'cloud_fraction2'   , cloud_frac2_10m    );
nc_varput ( ncf, 'cloud_fraction3'   , cloud_frac3_10m    );
nc_varput ( ncf, 'glitch_fraction'   , glitch_frac_10m    );
nc_varput ( ncf, 'cloud_fraction_threshold'   , cf_10m    ); % v1.0a
% NaNs seem to be handled internally, but not sure how universal this is

% global attributes regarding version, date, authorship
nc_attput ( ncf, nc_global, 'title', 'W-band 10-minute cloud fraction data for VOCALS 2008');
nc_attput ( ncf, nc_global, 'version', '1.0a' );
nc_attput ( ncf, nc_global, 'creation_date', datestr(datenum(clock),'yyyy-mm-dd') );
nc_attput ( ncf, nc_global, 'acknowledgement', 'Data are provided by Simon de Szoeke (Oregon State University), Chris Fairall (NOAA/ESRL/PSD), and Sandra Yuter (North Carolina State University) with funding from the NOAA Climate Program Office.');
nc_attput ( ncf, nc_global, 'author_contact_email', 'Simon de Szoeke <sdeszoek@coas.oregonstate.edu>');
% nc_attput ( ncf, nc_global, 'for_more_information', 'ftp://precip.meas.ncsu.edu/pub/vocals/VOCALS2008_WbandCloudFraction10min_1_0.readme.txt');

nc_addhist( ncf , 'v1.0a: threshold cloud fraction added')
