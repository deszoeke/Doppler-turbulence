% proc_wband_cloudtop.m
% 2009-04-08 :: VOCALS 2008 :: Simon de Szoeke
%
% Process the W-band radar data for cloud fraction stats.
%
% 2009-09-11 Update to include height correction for antenna zenith angle, when available. h=range*cos(zen_angle)
% 2011-02-22 Use absolute dBZ cloud fraction threshold

%cd /Users/sdeszoek/Data/cruises/VOCALS_2008/RHB/Scientific_analysis/programs/VOCALS2008_programs/wband
%addpath('~/Data/cruises/VOCALS_2008/RHB/Scientific_analysis/programs/VOCALS2008_programs/');
read_parameters;

% All radar data files
momentfile=dir([way_raw_data_wband '*MMCRMom.nc']);
year=2008;
base_time=NaN+zeros(length(momentfile),1);
base_time_offset=base_time;
ntimes=base_time;

% Memory strategy: 
% Cycle reflectivity data through 1 array.
% Write 1 file of 1-min cloud top for each MMCR file.
% Later concatenate files and average to 10-minute cloud top.

% get sizes to allocate Z array
for fi=1:length(momentfile)
    filename=[way_raw_data_wband momentfile(fi).name];
    A=nc_getvarinfo(filename,'time_offset');
    ntimes(fi)=A.Size(1);
end
height=nc_varget(filename,'Heights',[0 0],[1 -1]); % actually range
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
analog_noise_signal=-120.74; % dBmW peak diagnosed by Simon, includes proper offset already
digital_noise_signal=-115.4;  % dBmW
noise_margin=3.5; % dB
adhocthreshold=-43+dB_offset;
dBZnoiseavg=20*log10(h)+radar_const+analog_noise_signal;
dBZnoisefloor=20*log10(h)+radar_const+digital_noise_signal+noise_margin;
noisefloor=max(adhocthreshold,dBZnoisefloor);
% other additional thresholds for cloud and rain detection
%thr=[0 0.5 1 2 4 8 12]; %relative thresholds above noise
thr=-40:2:40; % absolute dBZ thresholds

% loop files
%starter=find(base_time_yday>=310,1,'first');
starter=460;
for fi=starter:length(momentfile) % fi>=2;
    %%prflname=[way_raw_data_wband momentfile(fi-1).name]; % previous file
    filename=[way_raw_data_wband momentfile(fi).name];
    
    % found no need to read the end of the previous file to get an integral minute
    % base_time in seconds since 1970-1-1 00:00:00
    base_time(fi)=nc_varget(filename,'base_time');
    base_time_offset(fi)=base_time(fi)-base_time(starter);
    % base_time_mld in matlab datenumber
    base_time_mld=double(base_time(fi))/86400 + datenum(1970,1,1,0,0,0);
    base_time_yday=base_time_mld-datenum(year,1,0,0,0,0);
    % time_offset in seconds
    time_offset=nc_varget(filename,'time_offset');
    % seconds since base_time(starter), for concatentation
    time_offset_cat=base_time_offset(fi)+time_offset;
    
    % reflectivity matrix
    %%Z(1:npre,:)=nc_varget_lapxm(prflname,'Reflectivity',[ntimes(fi-1)-npre 0],[npre length(h)]) ;
    %%Z(npre+1:npre+ntimes(fi),:)=nc_varget_lapxm(filename,'Reflectivity',[0 0],[-1 length(h)]);
    Z(1:ntimes(fi),:)=nc_varget_lapxm(filename,'Reflectivity',[0 0],[-1 length(h)])+dB_offset; %+dB_offset 2013-04-28 SPdeS (1.72 dB more liberal cloud detection was used for the 2011 paper)
    
    % Use a reflectivity criterion noisefloor defined above that selects clouds
    % and eliminates (almost) all noise returns. Further reduce 
    % nonmeteorological noise and clutter by requiring clouds to be
    % vertically homogeneous for 3 range gates.
    iscloud1=Z(1:ntimes(fi),:)-repmat(noisefloor,[ntimes(fi) 1])>0 & repmat(h,[ntimes(fi) 1])<2e3 & repmat(h,[ntimes(fi) 1])>400; % Z criterion for cloud
    iscloud2=iscloud1 & circshift(iscloud1,[0 1]); % also require cloud immediately below
    iscloud3=iscloud2 & circshift(iscloud1,[0 2]);
    iscloud=iscloud3;
    isglitch=iscloud1 & ~iscloud3; % 1- and 2-point clouds considered glitches
  
    if false
    % Compute cloud fraction using higher reflectivity thresholds (beyond
    % the noise floor)
    nz=length(h);
    nthr=length(thr);
    zthing=repmat(Z(1:ntimes(fi),:),[1 1 nthr]);
    tthing=shiftdim(repmat(thr',[1 ntimes(fi) nz 1]),1);
    nthing=repmat(noisefloor,[ntimes(fi) 1 nthr]);
    hthing=repmat(h,[ntimes(fi) 1 nthr]);
    ic=zthing>=tthing+nthing & hthing>400 & hthing<2100;
    end
    % Compute cloud fraction using absolute dBZ thresholds, and require Z>noise
    nz=length(h);
    nthr=length(thr);
    zthing=repmat(Z(1:ntimes(fi),:),[1 1 nthr]);
    tthing=shiftdim(repmat(thr',[1 ntimes(fi) nz 1]),1);
    nthing=repmat(noisefloor,[ntimes(fi) 1 nthr]);
    hthing=repmat(h,[ntimes(fi) 1 nthr]);
    ic=zthing>=tthing & zthing>=nthing & hthing>400 & hthing<2100;
    
    % Cloud fraction is determined for 60-second bins:
    time_offset_bin=floor(base_time(fi)/60)*60-base_time(fi) + (0:60:ceil(time_offset(end)/60)*60)';
    num_returns=zeros(length(time_offset_bin)-1,1);
    cldfrac1=NaN+zeros(length(time_offset_bin)-1,1);
    cldfrac2=cldfrac1;
    cldfrac3=cldfrac1;
    cf=repmat(cldfrac1,[1 length(thr)]);
    glitchfrac=cldfrac1;
    cloudtop_filealltime=zeros(size(time_offset));
    % loop 1 minute bins
    for i=1:length(time_offset_bin)-1
        in_bin=time_offset>=time_offset_bin(i) & time_offset<time_offset_bin(i+1);
        num_returns(i)=sum(in_bin);
        % simpler, more liberal cloud fraction counters
        cldfrac1(i)=sum(any(iscloud1(in_bin,:),2))/num_returns(i);
        cldfrac2(i)=sum(any(iscloud2(in_bin,:),2))/num_returns(i);
        cldfrac3(i)=sum(any(iscloud3(in_bin,:),2))/num_returns(i);
        glitchfrac(i)=sum(any(isglitch(in_bin,:),2))/num_returns(i);
        % higher thresholds
        cf(i,:)=sum(any(ic(in_bin,:,:),2))/num_returns(i);
    end % loop 1-min bins
    
    % write 1-minute cloud fraction and top height file for each MMCR file
    time_aggreg=base_time_offset(fi)+time_offset_bin(1:end-1); % seconds
    
%     fracfile=[way_proc_data_wband 'cloudfraction/' momentfile(fi).name(1:11) 'CloudFraction.txt']; % for simpler cloud counters
%     ff=fopen(fracfile,'w');
%     fprintf(ff,'%7.0f\t%6.4e\t%6.4e\t%6.4e\t%6.4e\n',[time_aggreg,cldfrac1,cldfrac2,cldfrac3,glitchfrac]');
%     fclose(ff);
    
    % for higher thresholds
    cfile=[way_proc_data_wband 'cloudfraction/' momentfile(fi).name(1:11) 'ClFrAbsThresh.txt'];
    ff=fopen(cfile,'w');
    fmt='%6.4e\t';
    fprintf(ff,['%7.0f\t' repmat(fmt,[1,length(thr)-1]) '%6.4e\n'],[time_aggreg,cf]');
    fclose(ff);
end % loop files

% if ~exist([way_proc_data_wband 'cloudfraction/CloudFraction2008310-336.txt'],'file')
%     filecat([way_proc_data_wband 'cloudfraction/2008*CloudFraction.txt'],...
%             [way_proc_data_wband 'cloudfraction/CloudFraction_1min_2008310-336.txt'])
% end
if ~exist([way_proc_data_wband 'cloudfraction/ClFrAbsThresh_1min_2008310-336.txt'],'file')

    filecat([way_proc_data_wband 'cloudfraction/2008*ClFrAbsThresh.txt'],...
            [way_proc_data_wband 'cloudfraction/ClFrAbsThresh_1min_2008310-336.txt'])
end

% load all the cloud top height data
A=load([way_proc_data_wband 'cloudfraction/CloudFraction_1min_2008310-336.txt']);
A(A==-999)=NaN; % missing value of -999 may be supplied (currently is not)
time=A(:,1);
cloud_frac1=A(:,2);
cloud_frac2=A(:,3);
cloud_frac3=A(:,4);
glitch_frac=A(:,5);

B=load([way_proc_data_wband 'cloudfraction/ClFrAbsThresh_1min_2008310-336.txt']);

if false
% PLOTS

%starter=1;
%base_time(starter)=1225916508;
time_yday=datenum(0,0,0,0,0,base_time(starter)+time)+datenum(1970,1,1,0,0,0)-datenum(2008,1,0,0,0,0);
bzz=0.004762*0.8;
subplot(2,1,1)
plot(time_yday,cloud_frac1+bzz*(0.5+rand(size(cloud_frac1))),'b.','markersize',1)
hold on
plot(time_yday,cloud_frac2+bzz*(0.5+rand(size(cloud_frac1))),'r.','markersize',1)
plot(time_yday,cloud_frac3+bzz*(0.5+rand(size(cloud_frac1))),'.','color',[0 .5 0],'markersize',1)
xlabel('2008 yearday')
ylabel('Cloud fraction')
title('VOCALS W-band radar cloud fraction')
axis([311 337 -0.1 1.1])
subplot(2,1,2)
plot(time_yday,cloud_frac1-cloud_frac2+bzz*(0.5+rand(size(cloud_frac1))),'b.','markersize',1)
hold on
plot(time_yday,cloud_frac2-cloud_frac3+bzz*(0.5+rand(size(cloud_frac1))),'r.','markersize',1)
%print('-dpng',[way_raw_images_wband 'CloudFraction_1min_wband.png'])

% compare the different cloud fraction estimates
bzzz=bzz*(0.5+rand(size(cloud_frac1)));
plot(nanmean(cloud_frac1),nanmean(cloud_frac2),'bo')
hold on
plot(nanmean(cloud_frac1),nanmean(cloud_frac3),'ro')
legend('2','3','Location','NorthWest')
plot(cloud_frac1+bzzz,cloud_frac2+bzz*(0.5+rand(size(cloud_frac1))),'.','markersize',1)
plot(cloud_frac1+bzzz,cloud_frac3+bzz*(0.5+rand(size(cloud_frac1))),'r.','markersize',1)
plot([0 1],[0 1],'k')
xlabel('single cloud')
ylabel('vertically contiguous cloud')
axis square
axis([0 1 0 1])
title '1-minute cloud fraction'
%print('-dpng',[way_proc_images_wband 'cloud_frac_1min_compare.png'])
end