% Load lookup table, hourly radar data, do corrections, perform Pinsky_retrieval_cloudrelative
% like call from Pinsky_main_cloudrelative.m but load all lookups for each hour (without ring memory)

%cd /Users/sdeszoek/Data/cruises/VOCALS_2008/RHB/Scientific_analysis/programs/VOCALS2008_programs/wband/Pinsky

% record-average lookup tables now inside 2nd huge loop
[Vgbin,Ustdbin,nbin]=deal(zeros(length(Zbin),length(hrel)));
[Vring,Uring,nring]=deal(zeros(length(Zbin),length(hrel),length(iw)));
wweight=repmat(shiftdim(weight,-1),[length(Zbin),length(hrel),1]); % temporal Gaussian weights

%% loop over files by day, hour
for iday=stday:enday; % back up 1 day to spin up lookup table?
    istatusday=find(iday==statusday); % index of motion status day (row)
    fprintf(1,'\n%3d ',iday)
    for idhr=0:23;
        istatushr=find(idhr==statushour); % index of motion status hour (col)
        fprintf(1,' %02d ',idhr)
        for iwi=1:25 % loop over hours to average in the window
            iwhr=mod(idhr+iw(iwi),24); % hours needed
            iwday=floor((iday*24+idhr+iw(iwi))/24);
            %% record-average lookup tables
            % composite time-local Z-hrel-bin lookup table from hourly tables
            fname=sprintf('/Users/sdeszoek/Data/cruises/VOCALS_2008/RHB/Scientific_analysis/programs/VOCALS2008_programs/wband/Pinsky/retrieval/Pinsky_lookup%04i_%03i_%02i_cloudrelative_p5dBZ.mat',...
                2008,iwday,iwhr); % the next hour needed for the window
            if exist(fname,'file')
                Tmp=load(fname);
                Tmp.Vgbin(isnan(Tmp.Vgbin))=0;
                Tmp.Ustdbin(isnan(Tmp.Ustdbin))=0;
                Tmp.nbin(isnan(Tmp.Vgbin))=0;
                Vring(:,:,iwi)=Tmp.Vgbin;
                Uring(:,:,iwi)=Tmp.Ustdbin;
                nring(:,:,iwi)=Tmp.nbin;
            else
                Vring(:,:,iwi)=0;
                Uring(:,:,iwi)=0;
                nring(:,:,iwi)=0;
            end
        end
        ssum=sum(nring.*wweight,3);
        Vgbin=sum(Vring.*nring.*wweight,3)./ssum;
        Ustdbin=sum(Uring.*nring.*wweight,3)./ssum;
        % robustness condition for lookup
        fac=sum(wweight(1,1,:),3)/size(wweight,3);
        Vgbin(ssum<16*fac)=NaN;
        Ustdbin(ssum<16*fac)=NaN;

        %% read moment and motion data, perform retrieval
        dddhh=sprintf('%03d%02d',iday,idhr);
        momentfile=dir([way_raw_data_wband '2008' dddhh '*MMCRMom.nc' ]);     
        kongfile=dir([way_raw_data_wband 'motion_adjT/' yyyy dddhh '*Kongsberg_adjT.txt']);
        crossfile=dir([way_raw_data_wband 'motion/' yyyy dddhh '*Crossbow.txt']);
        if ~isempty(momentfile) && ~isempty(kongfile)
            % read radar and motion data, synchronize, apply offsets, adjust Doppler vel.
            read_radar_motion
            
            % now perform the Pinsky retrieval
            [Vg,W,Ustd,a]=Pinsky_retrieval_cloudrelative(V,Z,Zbin,Vgbin,Ustdbin,rho,h,htop,hthick);
            U=V-Vg;
            Vgprime=U-W;
            % a=W./U;
            
            %% save Z, V(Doppler), Vg, Ustd, W
            time_offset=time_offset(1:countt);
            time_yday=time_yday_use(1:countt);
            fname=sprintf('/Users/sdeszoek/Data/cruises/VOCALS_2008/RHB/Scientific_analysis/programs/VOCALS2008_programs/wband/Pinsky/retrieval/Pinsky_wRetrieval%04i_%03i_%02i.mat',2008,iday,idhr);
            save(fname,'Z','V','Vg','Ustd','W','time_offset','time_yday','h')
            
        end % there is hourly data
    end % hr
end % day

return
