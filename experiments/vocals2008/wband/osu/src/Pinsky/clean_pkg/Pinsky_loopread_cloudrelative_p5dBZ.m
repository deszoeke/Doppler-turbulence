% Pinsky_loopread_cloudrelative_p5dBZ.m
% called by Pinsky_main_cloudrelative.m
% Load hourly radar data, do corrections, make lookup tables, save hourly lookup tables.
% calls Pinksy_makelookup_cloudrelative.m

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
        crossfile=dir([way_raw_data_wband 'motion/' yyyy dddhh '*Crossbow.txt']);
        if ~isempty(momentfile) && ~isempty(kongfile)
            % read radar and motion data, synchronize, apply offsets, adjust Doppler vel.
            read_radar_motion
            
            % make Pinsky retrieval step 1 lookup table (mean settling velocity component)
            [Vgbin,Ustdbin,nbin]=Pinsky_makelookup_cloudrelative(V,Z,Zbin,h,htop);          

            % save the hourly lookup tables
            save(sprintf('./retrieval/Pinsky_lookup%04i_%03i_%02i_cloudrelative_p5dBZ.mat',year,iday,idhr),...
                'h','Zbin','Vgbin','Ustdbin','nbin');

        end % there is hourly data
    end % hr
end % day

return
