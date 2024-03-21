% compare_top_sonde_wband
% VOCALS 2008 :: 2009-06-17 :: Simon de Szoeke
%
% Compare the cloud top from the W-band radar (1-minute resolution)
% and the inversion base from the soundings. They should coincide.

% get sounding
sonde=load('~/Data/cruises/VOCALS_2008/RHB/balloon/Processed/sonde_inversion.mat');

% get 1-minute w-band cloud top
A=load([way_proc_data_wband '/cloudheight/CloudHeight_1min_2008310-336.txt']);
A(A(:,3)<0,3)=NaN; % missing value of -999 for cloud top height may be supplied
time=A(:,1);
cloud_fraction=A(:,2);
cloud_top=A(:,3);
num_returns=A(:,4);
flag=A(:,5);
% flag values
% 0 no errors detected, cloud top computed.
% 1 cloud fraction=0, cloud top set to NaN.
% 2 cloud fraction>=0 but found no 3 cloud returns with consistent height,
%   possible glitch, provisional cloud top computed.

cloud_top(logical(flag))=NaN;
starter=1;
base_time(starter)=1225916508; % s since 1970-1-1 00:00
time_yday=datenum(0,0,0,0,0,base_time(starter)+time)+datenum(1970,1,1,0,0,0)-datenum(2008,1,0,0,0,0);

% find minute during which sonde was launched
minute=1/(60*24);
wband.cloudtop=zeros(length(sonde.yday),1);
wband.yday=wband.cloudtop;
wband.yday_mean=wband.cloudtop;
wband.yday_start=wband.cloudtop;
wband.yday_end=wband.cloudtop;
ideux=zeros(length(sonde.yday),1,'int16');
usesonde=logical(zeros(length(sonde.yday),1));
for isonde=1:length(sonde.yday)
    itest=find(time_yday<=sonde.yday(isonde),1,'last'); % index of wband corresp to sonde
    if ~isempty(itest)
        ideux(isonde)=itest;
        if sonde.yday(isonde)-time_yday(itest)>minute
            ideux(isonde)=0;
            usesonde(isonde)=0;
        else
            usesonde(isonde)=1;
            wband.cloudtop(isonde)=nanmean(cloud_top(itest-1:itest+8));
            wband.yday(isonde)=time_yday(itest);
            wband.yday_mean(isonde)=nanmean(time_yday(itest-1:itest+8));
            wband.yday_start(isonde)=time_yday(itest-1);
            wband.yday_end(isonde)=time_yday(itest+8);
        end
    else
        ideux(isonde)=0;
        usesonde(isonde)=0;
        continue
    end
    
end
iisonde=find(usesonde);

% plot
plot(time_yday,cloud_top/1e3,'.','markersize',3)
hold on
plot(sonde.yday,sonde.hinvbase/1e3,'r.')
plot(sonde.yday(iisonde),(sonde.hinvbase(iisonde)-wband.cloudtop(iisonde))/1e3,'k.-')
set(gca,'xlim',[315 337])
set(gca,'fontsize',14)
xlabel('2008 yearday')
ylabel('cloud top/inversion base height (km)')
legend('W-band','sonde','sonde-Wband')
print('-dpng',[way_proc_images_wband 'compare_cloudtop_sonde_wband.png']);
