% Pinsky_loopread_cloudrelative_1dBZ.m
% called by Pinsky_main_cloudrelative
% Load hourly radar data, do corrections, make lookup tables, save hourly lookup tables.
% do Pinksy_retrieval_cloudrelative

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
            Z=Z(1:countt,:); % truncate
            V=-V(1:countt,:); % positive upward
            
            % synchronize cloud top/base on radar time (yearday)
            htop=  interp1(synthtime,cloudtop  ,time_yday_use(1:countt),'nearest');
            hbase= interp1(synthtime,cloudbase ,time_yday_use(1:countt),'nearest');
            hthick=interp1(synthtime,cloudthick,time_yday_use(1:countt),'nearest');

            % make Pinsky retrieval step 1 lookup table (mean settling velocity component)
            [Vgbin,Ustdbin,nbin]=Pinsky_makelookup_cloudrelative(V,Z,Zbin,h,cloudtop4radar,cloudthick4radar);          
            % save the hourly lookup tables
            save(sprintf('./retrieval/Pinsky_lookup%04i_%03i_%02i_cloudrelative_1dBZ.mat',year,iday,idhr),...
                'h','Zbin','Vgbin','Ustdbin','nbin');
%             [Vgbin,Ustdbin,nbin]=Pinsky_makelookup_cloudrelative_thickness(V,Z,Zbin,h,htop,hthick,thickbin);
%             save(sprintf('./retrieval/Pinsky_lookup%04i_%03i_%02i_cloudrelative_thickness.mat',year,iday,idhr),...
%                 'h','Zbin','thickbin','Vgbin','Ustdbin','nbin'); 

        end % there is hourly data
    end % hr
end % day

return
