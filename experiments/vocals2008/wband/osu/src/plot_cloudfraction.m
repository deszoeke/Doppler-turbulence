% Plot 10-minute cloud fractions from ceilometer W-band radar.
% VOCALS 2008 :: 2009-10-07 :: Simon de Szoeke

read_parameters;
base_time=1225916508;

% 10-minute W-band cloud fractions
ncf=[way_proc_data_wband 'cloudfraction/VOCALS2008_WbandCloudFraction10min_1_0.nc'];
wcf.yday=nc_varget(ncf,'yday');
wcf.cf1=nc_varget(ncf,'cloud_fraction1');
wcf.cf2=nc_varget(ncf,'cloud_fraction2');
wcf.cf3=nc_varget(ncf,'cloud_fraction3');
wcf.gf=nc_varget(ncf,'glitch_fraction');

% 10-min ceilometer cloud fraction
A=load('/Users/sdeszoek/Data/cruises/VOCALS_2008/RHB/ceilometer/Processed/VOCALS2008ceilo_10min_a.txt');
ceil.yday=A(:,1);
ceil.cf=A(:,9);
ceil.height=A(:,13); % 85 percentile 1st cloud base height
ceil.cf(ceil.height>2000)=NaN; % filter above 2000m

% synchronize times to yday
yday=round(wcf.yday*60*24/10)/(60*24/10);
[iiw,locw]=ismember(round(wcf.yday*60*24/10),round(yday*60*24/10));
[iic,locc]=ismember(round(ceil.yday*60*24/10),round(yday*60*24/10));
% yday(locw(iiw)) approx equals wcf.yday(iiw)
wcf_yday=yday(locw(iiw));
wcf_1=wcf.cf1(iiw);
wcf_2=wcf.cf2(iiw);
wcf_3=wcf.cf3(iiw);
ceil_yday=yday(locc(iic));
ceil_cf=ceil.cf(iic);

[ii,loc]=ismember(wcf_yday,ceil_yday);
% ceil is shorter than wcf: wcf_yday(ii)=ceil_yday

% average to 1 hour
yday_h=ceil_yday(1:6:end); % ceil_yday(1) is 0300
wcf_1_h=ave_small2large(wcf_yday(ii),yday_h,wcf_1(ii));
wcf_2_h=ave_small2large(wcf_yday(ii),yday_h,wcf_2(ii));
wcf_3_h=ave_small2large(wcf_yday(ii),yday_h,wcf_3(ii));
ceil_cf_h=ave_small2large(ceil_yday,yday_h,ceil_cf);

% compare the different cloud fraction estimates (10 minute)
plot(nanmean(wcf_1),nanmean(wcf_2),'bo')
hold on
plot(nanmean(wcf_1),nanmean(wcf_3),'ro')
plot(nanmean(wcf_1),nanmean(ceil_cf),'mo')
legend('2','3','Location','NorthWest')
plot(wcf_1,wcf_2,'.','markersize',1)
plot(wcf_1,wcf_3,'r.','markersize',1)
plot(wcf_1(ii),ceil_cf,'m.','markersize',1)
plot([0 1],[0 1],'k')
xlabel('single cloud')
ylabel('vertically contiguous cloud')
axis square
axis([0 1 0 1])
title '10-minute cloud fraction'
% print('-dpng',[way_proc_images_wband 'cloud_frac_10min_compare.png'])

% compare the different cloud fraction estimates (1 minute)
plot(nanmean(wcf_1),nanmean(wcf_2),'bo')
hold on
plot(nanmean(wcf_1),nanmean(wcf_3),'ro')
plot(nanmean(wcf_1),nanmean(ceil_cf),'mo')
legend('2','3','ceil.','Location','NorthWest')
plot(wcf_1,wcf_2,'.','markersize',3)
plot(wcf_1,wcf_3,'r.','markersize',3)
plot(wcf_1(ii),ceil_cf,'mx','markersize',4)
plot([0 1],[0 1],'k')
xlabel('single W-band cloud')
ylabel('vertically contiguous/ceilometer cloud')
axis square
axis([0 1 0 1])
title '10-minute cloud fraction'
% print('-dpng',[way_proc_images_wband 'cloud_frac_10min_compare.png'])