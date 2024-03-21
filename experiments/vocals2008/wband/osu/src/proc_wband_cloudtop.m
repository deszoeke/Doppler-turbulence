% proc_wband_cloudtop.m
% 2009-04-08 :: VOCALS 2008 :: Simon de Szoeke
%
% Process the W-band radar data for cloud top.
%
% 2009-09-11 Update to include height correction for antenna zenith angle, when available. h=range*cos(zen_angle)

%cd /Users/sdeszoek/Data/cruises/VOCALS_2008/RHB/Scientific_analysis/programs/VOCALS2008_programs/wband
%addpath('~/Data/cruises/VOCALS_2008/RHB/Scientific_analysis/programs/VOCALS2008_programs/');
run('../read_parameters.m');

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
height=ncread(filename,'Heights',[1 1],[Inf 1]); % actually range
h=height(height<2e3)';                              %
range_gate=h(2)-h(1); % 24.9576 m
Z=NaN+zeros(max(ntimes),sum(height<2e3)); % preallocate, to be recycled

% After best recalibration, Ken Moran says to subtract 1.72 dB to all recorded values.
% reviever gain is 1.72 higher than previously estimated, which reduces dBZ
% and dBmW signals. This is an increase in sensitivity. Note recalibration noise
% source error of +-0.6 dB.
dB_offset=-1.72; % dB
% define a noise level for the reflectivity based on Ken Moran's manual
radar_const=19.613; % dB_
analog_noise_signal=-120.74; % dBmW noise peak diagnosed by Simon
digital_noise_signal=-115.4;  % dBmW
noise_margin=+3.5; % dB
adhocthreshold=-43+dB_offset;
dBZnoiseavg=20*log10(h)+radar_const+analog_noise_signal;
dBZnoisefloor=20*log10(h)+radar_const+digital_noise_signal+noise_margin;
noisefloor=max(adhocthreshold,dBZnoisefloor);

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
    Z(1:ntimes(fi),:)=nc_varget_lapxm(filename,'Reflectivity',[0 0],[-1 length(h)])+dB_offset;
    %+dB_offset 2013-04-28 SPdeS (1.72 dB more liberal cloud detection was used for the 2011 paper)
    
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
        H=repmat(h,[ntimes(fi) 1]); % no angle correction
    elseif length(kongfile)>1
        kongflag=2;
        fkongerr=fopen([way_raw_data_wband 'motion_adjT/kongsberg_log.txt'],'a+');
        fprintf(fkongerr,['Multiple files found: ' kongfile(:).name]);
        fclose(fkongerr);
        H=repmat(h,[ntimes(fi) 1]); % no angle correction
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
        %sintheta=sqrt(quad);
        %theta=asin(sintheta); % zenith angle (radians)
        costheta=sqrt(1-quad);
        if sum(isfinite(costheta))<max(3.5*60,length(costheta)/2)
            kongflag=3; % insufficient motion data
            H=repmat(h,[ntimes(fi) 1]); % no angle correction
            H(isfinite(costheta),:)=costheta(isfinite(costheta))*h; % correct where available
        else
            H=costheta*h; % antenna zenith angle cosine correction
        end
    end
    
    % Determine cloud top height for each return
    % vertical continuity BELOW is better than temporal here NOT DONE...
    if false % deglitch by requiring 50% cloud within 2.5 seconds (either side) of the return
        width=2.5; % s max half-width of window
        frac=0.30; % fraction of cloudy points to require in window
        [nt,nh]=size(iscloud2);
        iscloud=0*iscloud2;
        tic     % unf. slow...
        pickt=find(any(iscloud2,2));
        ismostlycloud=zeros(size(iscloud2));
        for ti=1:length(pickt)
            window=abs(time_offset-time_offset(pickt(ti)))<width;
            ismostlycloud(pickt(ti),:)=mean(iscloud2(window,:))>=frac;
        end
        iscloud=iscloud2 & ismostlycloud;
        iscloud(~ismostlycloud & iscloud2 & circshift(iscloud2,[0 1]))=1;
        toc
    end
    
%     %% encapsulate the rest of what's in this loop into a function
%     dt=60;
%     [flag,mean_top_bin,median_top_bin,val15_top_bin,val85_top_bin,cloud_fraction]=...
%         get_wband_cloudtop(Z(1:ntimes(fi),:),H,noisefloor,time_offset,dt,range_gate);
    
    % Use a reflectivity criterion noisefloor defined above that selects clouds
    % and eliminates (almost) all noise returns. Further reduce 
    % nonmeteorological noise and clutter by requiring clouds to be
    % vertically homogeneous for 3 range gates.
    iscloud1=Z(1:ntimes(fi),:)-repmat(noisefloor,[ntimes(fi) 1])>0 & repmat(h,[ntimes(fi) 1])<2e3 & repmat(h,[ntimes(fi) 1])>400; % Z criterion for cloud
    iscloud2=iscloud1 & circshift(iscloud1,[0 1]); % also require cloud immediately below
    iscloud3=iscloud2 & circshift(iscloud1,[0 2]);
    iscloud=iscloud3;
    isglitch=iscloud1 & ~iscloud3; % 1 and 2-point clouds considered glitches

    % cloud top height
    Cloudheight=H.*iscloud; % matrix of heights in cloud, 0 out of cloud, still has glitches
    soo=sort(Cloudheight(isfinite(Cloudheight)),'descend');
    highest=find(abs(soo(1:100)-soo(3:102))<range_gate/2,1,'first'); % highest cloud height for which there are at least 3 retrievals
    if ~isempty(highest)
        Cloudheight(Cloudheight>soo(highest))=0; % block too-high glitches
    else
        Cloudheight(:)=0; % cloudheight is all glitches
    end
    
    % Cloud top height is determined for 60-second bins:
    % The height for which there are at least 3 cloudy returns
    % in the minute is the highest allowed top. The cloud top for the
    % minute is taken as the mean of highest cloud tops at or below the
    % highest allowed top.
    time_offset_bin=floor(base_time(fi)/60)*60-base_time(fi) + (0:60:ceil(time_offset(end)/60)*60)';
    num_returns=zeros(length(time_offset_bin)-1,1);
    num_cloud_returns=zeros(length(time_offset_bin)-1,1);
    num_nonmet_returns=zeros(length(time_offset_bin)-1,1);
    flag=zeros(length(time_offset_bin)-1,1,'int8');
    mean_top_bin=NaN+zeros(length(time_offset_bin)-1,1);
    median_top_bin=mean_top_bin;
    std_top_bin=mean_top_bin;
    val15_top_bin=mean_top_bin;
    val85_top_bin=mean_top_bin;
    meanval_top_bin=mean_top_bin;
    cldfrac1=mean_top_bin;
    cldfrac2=mean_top_bin;
    cldfrac3=mean_top_bin;
    cloudtop_filealltime=zeros(size(time_offset));
    % loop 1 minute bins
    for i=1:length(time_offset_bin)-1
        in_bin=time_offset>=time_offset_bin(i) & time_offset<time_offset_bin(i+1);
        num_returns(i)=sum(in_bin);
        % is a cloud present?
        num_cloud_returns0=sum(any(iscloud(in_bin,:),2)); % cloud or nonmet
        % Separate cloud-or-nonmeteorological returns into cloud and nonmet by
        % insisting there are 3 returns at given height per minute for clouds.
        if num_cloud_returns0>0
            % cloud height
            C=Cloudheight(in_bin,:);
            soi=sort(C(:),'descend');
            highest=find(abs(soi(1:100)-soi(3:102))<range_gate/2 & soi(1:100)>0,1,'first'); % highest allowed: at least 3 at this level must be found
            if isempty(highest)
                flag(i)=2; % possible clouds, but no cloud with consistent height found
                num_cloud_returns(i)=0;
                num_nonmet_returns(i)=num_cloud_returns0; % assign all to nonmet
            else
                flag(i)=0; % no errors
                num_cloud_returns(i)=num_cloud_returns0-(highest-1);
                num_nonmet_returns(i)=highest-1;
                C(C>soi(highest))=0; % exclude tops above highest allowed                
            end
            % proceed to compute cloud top
            cloudtop_filealltime(in_bin)=max(C,[],2);
            valsort=sort(cloudtop_filealltime(in_bin & cloudtop_filealltime>0));
            lensort=sum(in_bin & cloudtop_filealltime>0);
            if lensort>=6 %%needs to be tested
                val15_top_bin(i)=valsort(round(0.15*lensort));
                val85_top_bin(i)=valsort(round(0.85*lensort));
                meanval_top_bin(i)=mean( valsort(round(0.15*lensort):round(0.85*lensort)) );
            else
                val15_top_bin(i)=NaN;
                val85_top_bin(i)=NaN;
                meanval_top_bin(i)=NaN;
            end
            median_top_bin(i)=median(cloudtop_filealltime(in_bin & cloudtop_filealltime>0));
            mean_top_bin(i)=mean(cloudtop_filealltime(in_bin & cloudtop_filealltime>0));
            std_top_bin(i)=std(cloudtop_filealltime(in_bin & cloudtop_filealltime>0));
        else % num_cloud_returns0==0
            flag(i)=1; % cloud fraction=0
            num_cloud_returns(i)=0;
            num_nonmet_returns(i)=0;
            cloudtop_filealltime(in_bin)=NaN;
            median_top_bin(i)=NaN;
            mean_top_bin(i)=NaN;
        end
        if num_cloud_returns<0
            disp('Whoa! Cloudy return number < 0')
            break
        end
                
        % simpler, more liberal cloud fraction counters
        cldfrac1(i)=sum(any(iscloud1(in_bin,:),2))/num_returns(i);
        cldfrac2(i)=sum(any(iscloud2(in_bin,:),2))/num_returns(i);
        cldfrac3(i)=sum(any(iscloud3(in_bin,:),2))/num_returns(i);
        glitchfrac(i)=sum(any(isglitch(in_bin,:),2))/num_returns(i);

    end % loop 1-min bins
    num_clear_returns=num_returns-num_cloud_returns-num_nonmet_returns;
    cloud_fraction=num_cloud_returns./(num_cloud_returns+num_clear_returns);
    
    % write 1-minute cloud fraction and top height file for each MMCR file
    time_aggreg=base_time_offset(fi)+time_offset_bin(1:end-1); % seconds
    %sprintf('%5.0f\t%6.4e\t%9.4f\t%3d\t%1d\n',[time_aggreg, cloud_fraction, mean_top_bin, double(num_returns), double(flag)]')
    
    heightfile=[way_proc_data_wband 'cloudheight/' momentfile(fi).name(1:11) 'CloudHeight.txt'];
    fh=fopen(heightfile,'w');
    fprintf(fh,['%7.0f\t%6.4e\t'...
                '%9.4f\t%9.4f\t'...
                '%9.4f\t%9.4f\t%9.4f\t'...
                '%9.4f\t'...
                '%3d\t%3d\t'...
                '%3d\t%3d\t'...
                '%1d\n'],...
        [time_aggreg, cloud_fraction,...
         max(-999,mean_top_bin), max(-999,std_top_bin),...
         max(-999,val15_top_bin), max(-999,median_top_bin), max(-999,val85_top_bin),...
         max(-999,meanval_top_bin),...
         double(num_returns), double(num_cloud_returns),...
         double(num_clear_returns), double(num_nonmet_returns),...
         double(flag)]');
    fclose(fh);
    
    fracfile=[way_proc_data_wband 'cloudfraction/' momentfile(fi).name(1:11) 'CloudFraction.txt']; % for simpler cloud counters
    ff=fopen(fracfile,'w');
    fprintf(ff,'%7.0f\t%6.4e\t%6.4e\t%6.4e\t%6.4e\n',[time_aggreg,cldfrac1,cldfrac2,cldfrac3,glitchfrac]');
    fclose(ff);
end % loop files

% flag values
% 0 no errors detected, cloud top computed.
% 1 cloud fraction=0, cloud top set to NaN.
% 2 cloud fraction>=0 but found no 3 cloud returns with consistent height,
%   possible glitch, provisional cloud top computed.

% NOTES
% Excluded 5 files
% 20083130129MMCRMom.nc
% 20083130200
% 20083131614
% 20083151432
% 20083151500
% in Raw/badcal because their reflectivities are 43 dBZ higher than normal.

% Cloud tops ~<=600 m at beginning and end of leg 2 (clutter: birds, bats, bugs,
% dust, haze? at coast). This is mostly solved by requiring cloudy returns
% to have 3 levels of vertically homogeneous Z above the cloud threshold.

% Cloud tops under Sc top are probably physical, but I don't know what they
% are. Could be rain not coincident with clouds above, rain attenuating
% radar, scud, or clutter.

% concatenate hourly 1-min data files if not done already
if ~exist([way_proc_data_wband 'cloudheight/CloudHeight2008310-336.txt'],'file')
    filecat([way_proc_data_wband 'cloudheight/2008*CloudHeight.txt'],...
            [way_proc_data_wband 'cloudheight/CloudHeight_1min_2008310-336.txt'])
end
if ~exist([way_proc_data_wband 'cloudfraction/CloudFraction2008310-336.txt'],'file')
    filecat([way_proc_data_wband 'cloudfraction/2008*CloudFraction.txt'],...
            [way_proc_data_wband 'cloudfraction/CloudFraction_1min_2008310-336.txt'])
end

% load all the cloud top height data
A=load([way_proc_data_wband 'cloudheight/CloudHeight_1min_2008310-336.txt']);
A(A==-999)=NaN; % missing value of -999 may be supplied
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
% nan out missing value=-999
cloud_top_val85(cloud_top_val85<0)=NaN;
cloud_top_val15(cloud_top_val15<0)=NaN;
cloud_top_valmean(cloud_top_valmean<0)=NaN;
% flag values
% 0 no errors detected, cloud top computed.
% 1 cloud fraction=0, cloud top set to NaN.
% 2 cloud fraction>=0 but found no 3 cloud returns with consistent height,
%   possible glitch, provisional cloud top computed.
% NaN out percentile-based cloud tops=0
lf=logical(flag);
cloud_top_val15(lf)=NaN;
cloud_top_val85(lf)=NaN;
cloud_top_valmean(lf)=NaN;
cloud_top_median(lf)=NaN;
cloud_top_std(lf)=NaN;

%starter=1;
%base_time(starter)=1225916508;
time_yday=datenum(0,0,0,0,0,base_time(starter)+time)+datenum(1970,1,1,0,0,0)-datenum(2008,1,0,0,0,0);
plot(time_yday,cloud_top_valmean/1e3,'r.-','markersize',3)
hold on
plot(time_yday,cloud_top_median/1e3,'b.','markersize',3)
set(gca,'fontsize',14)
xlabel('2008 yearday')
ylabel('Cloud top height (m)')
title('VOCALS W-band radar cloud top height')
axis([311 337 0 2])
print('-dpng',[way_raw_images_wband 'CloudHeight_1min_wband.png'])

