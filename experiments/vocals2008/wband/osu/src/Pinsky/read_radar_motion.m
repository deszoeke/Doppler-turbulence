countt=0; % number of data read so far
[Z,vel]=deal(NaN+zeros(12600,120));
for fi=1:size(momentfile,1) % handles multiple files
    filename=[way_raw_data_wband momentfile(fi,:).name];
    
    % no need to read the end of the previous file to get an integral minute
    % base_time in seconds since 1970-1-1 00:00:00
    base_time=ncread(filename,'base_time');
    % time_offset in seconds
    temp=ncread(filename,'time_offset');
    nt=length(temp);
    time_offset(countt+(1:nt),1)=temp;
    time(countt+(1:nt),1)=int32(time_offset(countt+(1:nt),1))+base_time;
    % base_time_mld in matlab datenumber
    base_time_mld=double(base_time)/86400 + datenum(1970,1,1,0,0,0);
    base_time_yday=base_time_mld-datenum(year,1,0,0,0,0);
    time_yday_use(countt+(1:nt),1)=base_time_yday+time_offset(countt+(1:nt),1)/86400;
    
    % reflectivity matrix
    Z(countt+(1:nt),:)=ncread(filename,'Reflectivity')'+dB_offset;
    %+dB_offset 2013-04-28 SPdeS (1.72 dB more liberal cloud detection was used for the 2011 paper);
    vel(countt+(1:nt),:)=ncread(filename,'MeanDopplerVelocity')';
    %width(countt+(1:nt),:)=ncread(filename,'SpectralWidth')';
    countt=countt+nt; % number of data read so far
end % multiple files

% Read Kongsberg motion compensation file
countk=0;
for ki=1:length(kongfile)
    kongfilename=[way_raw_data_wband 'motion_adjT/' kongfile(1).name]; % _adjT data starts on day 318
    [ktimetmp,kwtmp,pitchtmp,rolltmp]=read_kongsberg(kongfilename);
    nk=length(ktimetmp);
    kongtime(countk+(1:nk),1)=ktimetmp;
    kongw(countk+(1:nk),1)=kwtmp;
    pitch(countk+(1:nk),1)=pitchtmp;
    roll(countk+(1:nk),1)=rolltmp;
    countk=countk+nk;
end

% synchronize Kongsberg and radar time (yearday)
[temptime,ii,jj]=unique(kongtime(1:countk));
% do not interpolate over gaps more than 3 s
ibad=diff(temptime)>3/60/60/24;
ktilt=pitch(ii).^2+roll(ii).^2>10; % discard when antenna tilted--do not fill!
kongw(ii(ibad))=NaN;
kongw4radar=interp1(temptime,kongw(ii),time_yday_use(1:countt));

% backfill with Crossbow data if 10% Kong data missing and Crossbow status is good
if countk<0.9*36000 && Status(find(iday==statusday),idhr+1)==0
    countcross=0;
    for ci=1:length(crossfile)
        [crosstimetmp,crosswtmp]=read_crossbow([way_raw_data_wband 'motion/' crossfile(ci).name]);
        ncr=length(crosstimetmp);
        crosstime(countcross+(1:ncr),1)=crosstimetmp;
        crossw(countcross+(1:ncr),1)=crosswtmp;
        countcross=countcross+ncr;
    end
    ifill=isnan(kongw4radar);
    ii=isfinite(crossw(1:countcross)) & isfinite(crosstime(1:countcross));
    ibad=diff(crosstime(ii))>3/60/60/24;  % don't interpolate across >3s gap
    crossw(ii(ibad))=NaN;
    kongw4radar(ifill)=interp1(crosstime(ii),crossw(ii),time_yday_use(ifill));
end

% discard when antenna tilted
itilt=floor(unique(interp1(time_yday_use(1:countt),1:countt,temptime(ktilt))));
kongw4radar(itilt(isfinite(itilt)))=NaN;

% truncate data
Z=Z(1:countt,:);
vel=vel(1:countt,:);

% cutoff below noise floor
mask=Z>=repmat(noisefloor',[countt 1]);
Z(~mask)=NaN;
vel(~mask)=NaN;

% ship-motion-compensated vertical velocity +towards
V=vel(1:countt,:)+repmat(kongw4radar(1:countt),[1,length(h)]);

%width(~mask)=NaN;
% cutoff below 170 m range gate (clutter, near field contamination)
Z(:,h<170)=NaN;
V(:,h<170)=NaN;
%width(:,h<170)=NaN;
Z=Z(1:countt,:); % truncate
V=-V(1:countt,:); % positive upward

% synchronize cloud top/base on radar time (yearday)
htop=  interp1(time_yday_top,cloudtop  ,time_yday_use(1:countt),'nearest');
%             hbase= interp1(synthtime,cloudbase ,time_yday_use(1:countt),'nearest');
%             hthick=interp1(synthtime,cloudthick,time_yday_use(1:countt),'nearest');