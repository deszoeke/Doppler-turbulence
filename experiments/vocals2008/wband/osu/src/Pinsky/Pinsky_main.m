% Load hourly radar data, do corrections, feed to Pinsky_retrieval

%% set environment and constants
%cd /Users/sdeszoek/Data/cruises/VOCALS_2008/RHB/Scientific_analysis/programs/VOCALS2008_programs/wband/Pinsky
warning('off','MATLAB:interp1:NaNinY')
if exist('../../read_parameters.m','file')
    run('../../read_parameters')
else
    % W-band radar files
    way_raw_data_wband='/Users/sdeszoek/Data/cruises/VOCALS_2008/RHB/radar/wband/Raw/';
    way_proc_data_wband='/Users/sdeszoek/Data/cruises/VOCALS_2008/RHB/radar/wband/Processed/';
    way_raw_images_wband='/Users/sdeszoek/Data/cruises/VOCALS_2008/RHB/radar/wband/Raw_Images/';
    way_proc_images_wband='/Users/sdeszoek/Data/cruises/VOCALS_2008/RHB/radar/wband/Processed_Images/';
    year=2008;
end
yyyy=sprintf('%04d',year);
stday=318;
enday=336;
fs=3.5; %Hz

% parameters of Pinsky retrieval
Zbin=-40:0.25:10;
rho=0;
%rho=[0:64 repmat(64,1,55)]*.2/64; % approx correlation from Pinsky et al. Fig.5

momentfile=dir([way_raw_data_wband '2008*MMCRMom.nc' ]);
h=ncread([way_raw_data_wband momentfile(1).name],'Heights',[1 1],[Inf 1]);

%% get motion status, read from Ken and Sergio's spreadsheet
motion_status
% 0 OK
% 1 no data
% 2 bias (<1 degree)
% 3 noisy
% 4 motion adj. off or failed

%% After best recalibration, Ken Moran says to subtract 1.72 dB from all recorded values.
% receiver gain is 1.72 higher than previously estimated, which reduces dBZ
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

% loop length
nhr=24*(enday-stday+1);

%% loop over files by day, hour
ihr=0; % cumulative hour index
for iday=stday:enday;
    istatusday=find(iday==statusday); % index of motion status day (row)
    fprintf(1,'\n%3d ',iday)
    for idhr=0:23;
        istatushr=find(idhr==statushour); % index of motion status hour (col)
        fprintf(1,' %02d ',idhr)
        ihr=ihr+1;
        dddhh=sprintf('%03d%02d',iday,idhr);
        momentfile=dir([way_raw_data_wband '2008' dddhh '*MMCRMom.nc' ]);     
        kongfile=dir([way_raw_data_wband 'motion_adjT/' yyyy dddhh '*Kongsberg_adjT.txt']);
        if ~isempty(momentfile) && ~isempty(kongfile)
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
%                 time_mld=double(time)/86400 + datenum(1970,1,1,0,0,0);
%                 time_yday=time_mld-datenum(year,1,0,0,0,0);
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
                [ktimetmp,kwtmp]=read_kongsberg(kongfilename);
                nk=length(ktimetmp);
                kongtime(countk+(1:nk),1)=ktimetmp;
                kongw(countk+(1:nk),1)=kwtmp;
                countk=countk+nk;
            end
            
            % synchronize Kongsberg and radar time (yearday)
            [junk,ii,jj]=unique(kongtime);
            kongw4radar=interp1(kongtime(ii),kongw(ii),time_yday_use(1:countt));
                      
            % cutoff below noise floor
            mask=Z>=repmat(noisefloor',[length(Z) 1]);
            Z(~mask)=NaN;
            vel(~mask)=NaN;
            
            % ship-motion-compensated vertical velocity +towards
            V=vel(1:countt,:)+repmat(kongw4radar(1:countt),[1,length(h)]);

            %width(~mask)=NaN;
            % cutoff below 170 m range gate (clutter, near field contamination)
            Z(:,h<170)=NaN;
            V(:,h<170)=NaN;
            %width(:,h<170)=NaN;
            Z=Z(1:nt,:); % truncate
            V=-V(1:nt,:); % positive upward
            
            % call Pinsky retrieval (removes mean settling velocity component)
            [Vg,W,Vgbin,Ustdbin,Ustd,a]=Pinsky_retrieval(V,Z,Zbin);
            U=V-Vg;
            Vgprime=U-W;
            
            save(sprintf('./retrieval/Pinsky_lookup%i_%i_%i.mat',year,iday,idhr),...
                'h','Zbin','Vgbin','Ustdbin'); % just save the lookup tables
%             % plots
%             Pinsky_plot
        end % there is hourly data
    end % hr
end % day