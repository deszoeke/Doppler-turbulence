% Combine 10-minute cloud boundary data from W-band, ceilometer, and soundings
% write version 1.0
% VOCALS 2008 :: 2009-12-30 :: Simon de Szoeke

read_parameters;
ver=1.0;

%load([way_proc_data_wband '1min_stat/Z_1min.mat']);
%load([way_proc_data_wband '1min_stat/w_1min.mat']);

% load 10-min W-band cloud top height (Simon de Szoeke, Oregon State U.)
chf='/Users/sdeszoek/Data/cruises/VOCALS_2008/RHB/radar/wband/Processed/cloudheight/VOCALS2008_WbandCloudHeight10min_1_0.nc';
cth10.yday=nc_varget(chf,'yday');
cth10.height=nc_varget(chf,'cloudtop');
cth10.frac=nc_varget(chf,'cloudfrac');
cth10.numreturns=nc_varget(chf,'numreturns');
cth10.numcloud=nc_varget(chf,'numcloud');
cth10.numclear=nc_varget(chf,'numclear');
cth10.numnonmet=nc_varget(chf,'numnonmet');

% load 10-min CL31 cloud base height and cloud fraction
A=load('/Users/sdeszoek/Data/cruises/VOCALS_2008/RHB/ceilometer/Processed/VOCALS2008ceilo_10min_a.txt');
clb10.yday=A(:,1);
clb10.nobsall=A(:,2);
clb10.nobsclear=A(:,3);
clb10.nobscloud=A(:,4)+A(:,5);
clb10.nobsobsc=A(:,6)+A(:,7);
clb10.clearfrac=A(:,8);
clb10.cloudfrac=A(:,9);
clb10.height=A(:,13); % 85 percentile 1st cloud base height
clb10.height(clb10.height>2000)=NaN; % filter above 2000m
% load cb, ceilometer-bacscatter threshold cloud fraction (10-min)
load([way_proc_data_ceilo 'cloudfrac_backscatterthreshold.mat']);
%cb.yday == wcf.yday

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

% load inversion heights
invfile=[way_proc_data_balloon 'VOCALS2008_sonde_inversion.nc'];
sonde_yday=nc_varget(invfile,'yday');
sonde_lat=nc_varget(invfile,'lat');
sonde_lon=nc_varget(invfile,'lon');
%hinvtop=nc_varget(invfile,'height_inv_top');
hinvbase=nc_varget(invfile,'height_inv_base');
%tinvtop=nc_varget(invfile,'temp_inv_top');
%tinvbase=nc_varget(invfile,'temp_inv_base');
%thinvtop=nc_varget(invfile,'theta_inv_top');
%thinvbase=nc_varget(invfile,'theta_inv_base');
%dthdzinv=(thinvtop-thinvbase)./(hinvtop-hinvbase);

% synchronize times (10 min)
time10yday=( floor(min(sonde_yday)*24*60/10)*10/(24*60):10/(24*60):...
             max([clb10.yday;sonde_yday])+5/(24*60) )'; % even 10-min
n10m=length(time10yday);
% cloudfrac10a=nanmean(reshape(cloudfrac1,[10 n10m]))';

% synchronize LWP
[ii,loc]=ismember(round(lwp.yday*60*24/10),round(time10yday*60*24/10));
% time10yday(loc(ii)) approx equals lwp.yday(ii)
lwp10=NaN+time10yday;
adlwp10=lwp10;
lwp10(loc(ii))=lwp.lwp(ii);
adlwp10(loc(ii))=lwp.adlwp(ii);

% synchronize W-band cloud top height cth10 to cloudtop10(time10yday)
[ii,loc]=ismember(round(cth10.yday*60*24/10),round(time10yday*60*24/10));
cloudtop10=NaN+time10yday;
numreturns10=cloudtop10;
numclear10=cloudtop10;
numcloud10=cloudtop10;
numnonmet10=cloudtop10;
cloudtop10(loc(ii))=cth10.height(ii);
% cloudfrac10=NaN+time10yday;
% cloudfrac10(loc(ii))=cth10.frac(ii);
% cloudfracth10=NaN+zeros(max(loc),length(wcf.thr));
% cloudfracth10(loc(ii),:)=wcf.cft(ii,:);
numreturns10(loc(ii))=cth10.numreturns(ii);
numclear10(loc(ii))=cth10.numclear(ii);
numcloud10(loc(ii))=cth10.numcloud(ii);
numnonmet10(loc(ii))=cth10.numnonmet(ii);

%synchronize ceilo cloud base height clb10 to cloudbase10(time10yday)
[ii,loc]=ismember(round(clb10.yday*60*24/10),round(time10yday*60*24/10));
cloudbase10=NaN+time10yday;
nobsall10=NaN+time10yday;
nobscloud10=NaN+time10yday;
nobsclear10=NaN+time10yday;
nobsobsc10=NaN+time10yday;
cloudbase10(loc(ii))=clb10.height(ii);
nobsall10(loc(ii))=clb10.nobsall(ii);
nobscloud10(loc(ii))=clb10.nobscloud(ii);
nobsclear10(loc(ii))=clb10.nobsclear(ii);
nobsobsc10(loc(ii))=clb10.nobsobsc(ii);
% also sync ceilometer cloud fraction ceilfrac10(time10yday)
ceilfrac10=NaN+time10yday;
ceilfrac10(loc(ii))=clb10.cloudfrac(ii);

% interpolate lat and lon
lat10=interp1(flux.yday,flux.lat,time10yday);
lon10=interp1(flux.yday,flux.lon,time10yday);

% interpolate cloud top from soundings
sondecloudtop10=interp1(sonde_yday,hinvbase,time10yday);
sondeflag=zeros(size(sondecloudtop10));
sondeflag(isfinite(sondecloudtop10))=2;
[lognear,ind]=ismember(round(time10yday*24*60/10),round(sonde_yday*24*60/10));
nearst=find(lognear);
% for i=1:length(sonde_yday)
%     nearst(i,1)=find(sonde_yday(i)<time10yday(2:end) & sonde_yday(i)>=time10yday(1:end-1));
% end
sondecloudtop10(nearst)=hinvbase;
sondeflag(nearst)=1;
gap=find(diff(sonde_yday)>.25); % do not interpolate across gaps of more than 6 hours
for i=1:length(gap)
    ii=time10yday>sonde_yday(gap(i)) & time10yday<sonde_yday(gap(i)+1)-10/(24*60);
    sondecloudtop10(ii)=NaN;
    sondeflag(ii)=0;
end
% sondeflag
% 0 not available
% 1 collocated with sounding
% 2 interpolated from soundings

% QC to zap
% cloudbase above any wband cloud top or 100 m above interpolated 
cloudbase10(cloudbase10>cloudtop10)=NaN;
% wband cloudtop below cloud base AND below interpolated sonde cloud top
cloudtop10(cloudbase10>cloudtop10 & cloudtop10<sondecloudtop10)=NaN;
% cloud base 100 m above interpolated sonde cloud top
cloudbase10(cloudbase10-sondecloudtop10>100)=NaN;
% ashore between legs 1 and 2
cloudtop10(time10yday>309 & time10yday<315)=NaN;
cloudbase10(time10yday>309 & time10yday<315)=NaN;

% Write 10 minute W-band CLOUD TOP height in a netcdf file.
% write sounding INVERSION BASE height.
% write a synthesis CLOUD TOP that uses 10-min W-band where approporiate and ~4 hourly soundings where appropriate.
% write 10 minute ceilometer CLOUD BASE.
ncpath='/Users/sdeszoek/Data/cruises/VOCALS_2008/RHB/radar/wband/Processed/cloudheight/';
% note netcdflib does not accept ~/ as a path
ncf=[ncpath 'VOCALS2008CloudBoundaries10min_v' sprintf('%3.1f',ver) '.nc'];
nc_create_empty( ncf );
nc_add_dimension( ncf, 'time', 0 ) % record dimension

varstruct.Name = 'yday';
varstruct.Nctype = nc_float;
varstruct.Dimension = {'time'}; % { cell } needed
nc_addvar ( ncf , varstruct )
% more floats...
varstruct.Name = 'lat'; nc_addvar ( ncf , varstruct );
varstruct.Name = 'lon'; nc_addvar ( ncf, varstruct );
varstruct.Name = 'cloudtop'; nc_addvar ( ncf , varstruct );
varstruct.Nctype = nc_int; % some integers...
varstruct.Name = 'numwband'; nc_addvar ( ncf , varstruct );
varstruct.Name = 'numwbandcloud'; nc_addvar ( ncf , varstruct );
varstruct.Name = 'numwbandclear'; nc_addvar ( ncf , varstruct );
varstruct.Name = 'numwbandnonmet'; nc_addvar ( ncf , varstruct );
varstruct.Nctype = nc_float;
varstruct.Name = 'cloudbase'; nc_addvar ( ncf , varstruct);
varstruct.Nctype=nc_int;
varstruct.Name = 'numceil'; nc_addvar ( ncf , varstruct);
varstruct.Name = 'numceilcloud'; nc_addvar ( ncf , varstruct);
varstruct.Name = 'numceilclear'; nc_addvar ( ncf , varstruct);
varstruct.Name = 'numceilobscured'; nc_addvar ( ncf , varstruct);
varstruct.Nctype = nc_float;
varstruct.Name = 'sondecloudtop'; nc_addvar ( ncf , varstruct );
varstruct.Nctype=nc_int;
varstruct.Name = 'sondeinterpflag'; nc_addvar ( ncf , varstruct);

% variable attributes
nc_attput ( ncf , 'yday', 'long_name', 'yearday since 2007-December-31 00:00 (1.0 for 2008 January 1 00:00)' );
nc_attput ( ncf , 'yday', 'units', 'day' );
nc_attput ( ncf , 'yday', 'description', 'beginning of time interval');
nc_attput ( ncf , 'lat', 'units', 'degrees north' );
nc_attput ( ncf , 'lon', 'units', 'degrees east' );

nc_attput ( ncf , 'cloudtop', 'long_name', '10-minute mean W-band radar mid-range-gate height of highest cloud');
nc_attput ( ncf , 'cloudtop', 'units', 'm');
nc_attput ( ncf , 'numwband', 'long_name', 'total number of W-band returns in 10-minute time interval');
nc_attput ( ncf , 'numwband', 'units', 'unitless number');
nc_attput ( ncf , 'numwband', 'formula', 'numwband = numwbandcloud + numwbandclear + numwbandnonmet');
nc_attput ( ncf , 'numwbandcloud', 'long_name', 'number of cloudy returns in time interval');
nc_attput ( ncf , 'numwbandcloud', 'units', 'unitless number');
nc_attput ( ncf , 'numwbandcloud', 'note', 'W-band and ceilometer have different sensitivity to clouds')
nc_attput ( ncf , 'numwbandclear', 'long_name', 'number of clear (to W-band) returns in time interval');
nc_attput ( ncf , 'numwbandclear', 'units', 'unitless number');
nc_attput ( ncf , 'numwbandnonmet', 'long_name', 'number of nonmeteorological returns in time interval');
nc_attput ( ncf , 'numwbandnonmet', 'units', 'unitless number');
nc_attput ( ncf , 'numwbandnonmet', 'description', 'questionable cloud of limited vertical extent, noise, or clutter suspected')
nc_attput ( ncf , 'numwbandnonmet', 'long_name', 'number of nonmeteorological returns in time interval');
nc_attput ( ncf , 'numwbandnonmet', 'units', 'unitless number');
nc_attput ( ncf , 'numwbandnonmet', 'description', 'questionable cloud of limited vertical extent, noise, or clutter suspected')
nc_attput ( ncf , 'cloudbase', 'long_name', '10-minute 85th percentile (highest) ceilometer cloud base');
nc_attput ( ncf , 'cloudbase', 'units', 'm');
nc_attput ( ncf , 'numceil', 'long_name', 'total number of ceilometer returns in 10-minute time interval');
nc_attput ( ncf , 'numceil', 'units', 'unitless number');
nc_attput ( ncf , 'numceil', 'formula', 'numceil = numceilcloud + numceilclear + numceilobscure');
nc_attput ( ncf , 'numceilcloud', 'long_name', 'number of cloudy returns in time interval');
nc_attput ( ncf , 'numceilcloud', 'units', 'unitless number');
nc_attput ( ncf , 'numceilcloud', 'note', 'W-band and ceilometer have different sensitivity to clouds')
nc_attput ( ncf , 'numceilclear', 'long_name', 'number of clear (to W-band) returns in time interval');
nc_attput ( ncf , 'numceilclear', 'units', 'unitless number');
nc_attput ( ncf , 'numceilobscured', 'long_name', 'number of obscured returns in time interval');
nc_attput ( ncf , 'numceilobscured', 'units', 'unitless number');
nc_attput ( ncf , 'numceilobscured', 'description', 'ceilometer identifies obscured or partially obscured sky')
nc_attput ( ncf , 'sondecloudtop', 'long_name','radiosonde (~4hour) cloud top interpolated to 10-minutes');
nc_attput ( ncf , 'sondecloudtop', 'units', 'm');
nc_attput ( ncf , 'sondeinterpflag', 'long_name','interpolation from radiosonde status flag');
nc_attput ( ncf , 'sondeinterpflag', 'description', '0=not available; 1=collocated with sounding; 2=interpolated from sounding');

% global attributes regarding version, date, authorship
nc_attput ( ncf, nc_global, 'title', '10-minute cloud top height and base height data for VOCALS 2008');
nc_attput ( ncf, nc_global, 'version', sprintf('%3.1f',ver) );
nc_attput ( ncf, nc_global, 'creation_date', datestr(datenum(clock),'yyyy-mm-dd') );
nc_attput ( ncf, nc_global, 'acknowledgement', 'Data are provided by Simon de Szoeke (Oregon State University), Chris Fairall and Daniel Wolfe (NOAA/ESRL/PSD), and Sandra Yuter (North Carolina State University) with funding from the NOAA Climate Program Office.');
nc_attput ( ncf, nc_global, 'author_contact_email', 'Simon de Szoeke <sdeszoek@coas.oregonstate.edu>');
nc_attput ( ncf, nc_global, 'for_more_information', ['ftp://ftp.coas.oregonstate.edu/dist/sdeszoek/vocals/VOCALS2008CloudBoundaries10min_v' sprintf('%3.1f',ver) '.readme.txt']);

% write the data
nc_varput ( ncf, 'yday'       , time10yday  );
nc_varput ( ncf, 'lat'   , lat10       );
nc_varput ( ncf, 'lon'  , lon10       );
nc_varput ( ncf, 'cloudtop'   , cloudtop10   );
nc_varput ( ncf, 'numwband' , numreturns10 );
nc_varput ( ncf, 'numwbandcloud'  , numcloud10   );
nc_varput ( ncf, 'numwbandclear'  , numclear10   );
nc_varput ( ncf, 'numwbandnonmet' , numnonmet10  );
nc_varput ( ncf, 'cloudbase'      , cloudbase10 );
nc_varput ( ncf, 'numceil'       , nobsall10);
nc_varput ( ncf, 'numceilcloud'  , nobscloud10);
nc_varput ( ncf, 'numceilclear'  , nobsclear10);
nc_varput ( ncf, 'numceilobscured' , nobsobsc10);
nc_varput ( ncf, 'sondecloudtop' , sondecloudtop10);
nc_varput ( ncf, 'sondeinterpflag' , sondeflag);
% NaNs seem to be handled internally, but not sure how universal this is

%%%%
if false
rnd=-5+10*rand(size(cloudbase10));
plot(cloudbase10+rnd,cloudtop10,'-','linewidth',.1)
hold on
plot(cloudbase10+rnd,cloudtop10,'.','markersize',3)
plot([600 2000],[600 2000],'k')
plot([600 1300],[1300 2000],'k:')
axis([600 2000 600 2000])
axis square
set(gca,'fontsize',18)
xlabel('cloud base height')
ylabel('cloud top height')
text(640,1900,'rain & Cu under Sc')
print('-dpng',[way_proc_images_wband 'cloud_top_vs_base.png'])
print('-depsc2',[way_proc_images_wband 'cloud_top_vs_base.eps'])

end