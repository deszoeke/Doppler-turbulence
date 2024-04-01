% proc_wband_cloudtop_10min
% VOCALS 2008 :: 2009-06-19 :: Simon de Szoeke
%
% Average cloud top height from W-band radar (1-minute resolution)
% to 10 minute resolution. Requires output from proc_wband_cloudtop.m
% to be available first.

% preamble
run('../read_parameters')

% load all the cloud fraction and top height data
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

% AVERAGE to 10-minute (600 s) interval
sec1m=60;
sec10min=60*10;
% varaibles are named: name_units_resolution: sec=seconds, 1m=1 minute, 10m=10 minutes
time_sec_1m=base_time(starter)+time;
time0_10m_1m=floor((base_time(starter)+time)/sec10min);
time0_10m_10m=unique(time0_10m_1m);
time0_sec_10m=unique(time0_10m_1m)*sec10min; %starting time in seconds of 10m interval
n10m=length(unique(time0_10m_10m));
cloudtop_10m=NaN+zeros(n10m,1);
cloudfrac_10m=NaN+zeros(n10m,1);
numcloud_10m=NaN+zeros(n10m,1);
numclear_10m=NaN+zeros(n10m,1);
numnonmet_10m=NaN+zeros(n10m,1);
numreturns_10m=NaN+zeros(n10m,1);
for ti=1:n10m
    iitime=time0_10m_1m==time0_10m_10m(ti);
    ii=(iitime & isfinite(cloud_top_mean));
    temptop=sort(cloud_top_mean(ii)); % meters
    cloudtop_10m(ti)=mean(temptop); % straight mean
    cloudfrac_10m(ti)=mean(cloud_fraction(iitime)); % weights each minute the same
    numcloud_10m(ti)=sum(num_cloud_returns(ii)); % number of cloudy returns
    numclear_10m(ti)=sum(num_clear_returns(ii));
    numnonmet_10m(ti)=sum(num_nonmet_returns(ii));
    numreturns_10m(ti)=sum(num_returns(iitime)); % total number of returns used to compute the cloud fraction
end
time0_yday_10m=datenum(0,0,0,0,0,time0_sec_10m)+datenum(1970,1,1,0,0,0)-datenum(2008,1,0,0,0,0);

% plot to check
plot(time0_yday_10m,cloudtop_10m,'r.')
hold on
plot(time_yday,cloud_top_mean,'.','markersize',3)

% save the cloud top in a text file
heightfile=[way_proc_data_wband 'cloudheight/VOCALS2008_WbandCloudHeight10min_1_0.txt'];
fh=fopen(heightfile,'w');
%fprintf(fh,'%12.8f\t%10.3f\t%10.5f\t%4d\t%4d\t%4d\t%4d\n',...
fprintf(fh,'%12.8f\t%10.3f\t%10.5f\t%4d\t%4d\t%4d\t%4d\n',...
    [time0_yday_10m, max(-999,cloudtop_10m), max(-999,cloudfrac_10m),...
    numreturns_10m numcloud_10m numclear_10m numnonmet_10m]');
fclose(fh);

% to read the data
%heightfile=[way_proc_data_wband 'cloudheight/CloudHeight_10min.txt'];
missing_value=-999;
fh=fopen(heightfile);
B=fscanf(fh,'%f',[7 3234])';
B(B==missing_value)=NaN;
yday=B(:,1);
cloudtop=B(:,2);
cloudfrac=B(:,3);
numreturns=B(:,4);
numcloud=B(:,5);
numclear=B(:,6);
numnonmet=B(:,7);
fclose(fh);
% time0_yday_10m=yday;
% cloudtop_10m=cloudtop;
% cloudfrac_10m=cloudfrac;
% numreturns_10m=numreturns;
% numcloud_10m=numcloud;
% numclear_10m=numclear;
% numnonmet_10m=numnonmet;

% save the cloud top in a netcdf file
ncpath='/Users/sdeszoek/Data/cruises/VOCALS_2008/RHB/radar/wband/Processed/cloudheight/';
% note netcdflib does not accept ~/ as a path
ncf=[ncpath 'VOCALS2008_WbandCloudHeight10min_1_0.nc'];
nc_create_empty( ncf );
nc_add_dimension( ncf, 'time', 0 ) % record dimension

varstruct.Name = 'yday';
varstruct.Nctype = nc_float;
varstruct.Dimension = {'time'}; % { cell } needed
nc_addvar ( ncf , varstruct )
% more floats...
varstruct.Name = 'cloudtop'; nc_addvar ( ncf , varstruct );
varstruct.Name = 'cloudfrac'; nc_addvar ( ncf , varstruct );
varstruct.Nctype = nc_int; % some integers...
varstruct.Name = 'numreturns'; nc_addvar ( ncf , varstruct );
varstruct.Name = 'numcloud'; nc_addvar ( ncf , varstruct );
varstruct.Name = 'numclear'; nc_addvar ( ncf , varstruct );
varstruct.Name = 'numnonmet'; nc_addvar ( ncf , varstruct );

% variable attributes
nc_attput ( ncf , 'yday', 'long_name', 'yearday since 2007-December-31 00:00 (1.0 for 2008 January 1 00:00)' );
nc_attput ( ncf , 'yday', 'units', 'day' );
nc_attput ( ncf , 'yday', 'description', 'beginning of time interval');
nc_attput ( ncf , 'cloudtop', 'long_name', 'mean mid-range-gate height of highest cloud');
nc_attput ( ncf , 'cloudtop', 'units', 'm');
%nc_attput ( ncf , 'cloudtop', 'missing_value', missing_value)

nc_attput ( ncf , 'cloudfrac', 'long_name', 'cloud fraction');
nc_attput ( ncf , 'cloudfrac', 'units', 'unitless fraction');
nc_attput ( ncf , 'cloudfrac', 'description', 'fraction of returns in time interval with vertically coherent Z above noise threshold');
%nc_attput ( ncf , 'cloudfrac', 'missing_value', missing_value)

nc_attput ( ncf , 'numreturns', 'long_name', 'total number of returns in time interval');
nc_attput ( ncf , 'numreturns', 'units', 'unitless number');
nc_attput ( ncf , 'numreturns', 'formula', 'numreturns = numcloud + numclear + numnonmet');

nc_attput ( ncf , 'numcloud', 'long_name', 'number of cloudy returns in time interval');
nc_attput ( ncf , 'numcloud', 'units', 'unitless number');

nc_attput ( ncf , 'numclear', 'long_name', 'number of clear returns in time interval');
nc_attput ( ncf , 'numclear', 'units', 'unitless number');

nc_attput ( ncf , 'numnonmet', 'long_name', 'number of nonmeteorological returns in time interval');
nc_attput ( ncf , 'numnonmet', 'units', 'unitless number');
nc_attput ( ncf , 'numnonmet', 'description', 'questionable cloud of limited vertical extent, noise, or clutter suspected')

% write the data
nc_varput ( ncf, 'yday'       , yday       );
nc_varput ( ncf, 'cloudtop'   , cloudtop   );
nc_varput ( ncf, 'cloudfrac'  , cloudfrac  );
nc_varput ( ncf, 'numreturns' , numreturns );
nc_varput ( ncf, 'numcloud'   , numcloud   );
nc_varput ( ncf, 'numclear'   , numclear   );
nc_varput ( ncf, 'numnonmet'  , numnonmet  );
% NaNs seem to be handled internally, but not sure how universal this is

% global attributes regarding version, date, authorship
nc_attput ( ncf, nc_global, 'title', 'W-band 10-minute cloud fraction and cloud top height data for VOCALS 2008');
nc_attput ( ncf, nc_global, 'version', '1.0' );
nc_attput ( ncf, nc_global, 'creation_date', datestr(datenum(clock),'yyyy-mm-dd') );
nc_attput ( ncf, nc_global, 'acknowledgement', 'Data are provided by Simon de Szoeke (Oregon State University), Chris Fairall (NOAA/ESRL/PSD), and Sandra Yuter (North Carolina State University) with funding from the NOAA Climate Program Office.');
nc_attput ( ncf, nc_global, 'author_contact_email', 'Simon de Szoeke <sdeszoek@coas.oregonstate.edu>');
nc_attput ( ncf, nc_global, 'for_more_information', 'ftp://ftp.coas.oregonstate.edu/dist/sdeszoek/vocals/VOCALS2008_WbandCloudHeight10min_1_0.readme.txt');
