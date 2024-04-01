% compile_wband_1min_stat.m
% VOCALS 2008 :: 2009-07-07 :: Simon de Szoeke
%
% Compile 5 vertically-resolved statistics of Z and w:
% mean, min, max, variance, and skewness
%
% Compile vertically-resolved histograms of reflectivity, Doppler velocity w,
% and Doppler width (v.6).
%
% Input serial unformatted binary files from proc_wband_1min_stat_v6.m.
% Compile cloud top into a netcdf file with proc_wband_cloudtop_10min.m.
%
% Dimensions:
% time( )                               serial time (s), record dimension
% height(75)    134:1982                center of range gates (m)
% Zbins(38)     [-Inf -40:2:30 Inf]     reflectivity histogram bin edges (dBZ)
% wbins(53)     [-Inf -5:0.2:5 Inf]     velocity histogram bins edges (m/s)
% Dwbins(73)    [-Inf  0:0.1:7 Inf];    Doppler vel width hist bin edges (m/s?)

read_parameters;
load([way_proc_data_wband '1min_stat/Dw_1min.mat'], 'Dw','time_aggreg','time_yday');
% now fix Doppler width cfad

% % concatenate 1-min data files if not done already
% if ~exist([way_proc_data_wband '1min_stat/VOCALS2008_1min_stat.dat'],'file')
%     filecat([way_proc_data_wband '1min_stat/2008*_1min_stat.dat'],...
%             [way_proc_data_wband '1min_stat/VOCALS2008_1min_stat.dat'])
% end
if ~exist([way_proc_data_wband '1min_stat/onlywidth/VOCALS2008_cfad.dat'],'file')
    filecat([way_proc_data_wband '1min_stat/onlywidth/2008*_cfad.dat'],...
            [way_proc_data_wband '1min_stat/onlywidth/VOCALS2008_cfad.dat'])
end

% Define dimensions
% dBZ, Doppler vel, and width bin edges
% Zbins=[-inf -40:2:30 inf];
% wbins=[-inf -5:0.2:5 inf];
Dwbins=[-inf 0:0.1:7 inf];
% vertical coordinate (center of ranges, m)
h=load([way_proc_data_wband '1min_stat/range.txt']);
nheight=length(h); % 75
nday=25;
nmin=60;
nhour=24;
% nZbin=length(Zbins);
% nwbin=length(wbins);
nDwbin=length(Dwbins);

% %mean(time[ ],height[75])
% time_aggreg=zeros(nmin*nhour*nday,1);
% kongflag=time_aggreg; % V6 2010-07-13
% zen_angle.mean=NaN+zeros(nmin*nhour*nday,1);
% zen_angle.std=zen_angle.mean;
% zen_angle.cos=zen_angle.mean;
% zen_angle.sin=zen_angle.mean;
% %A=[Z.mean Z.min Z.max Z.var Z.skew w.mean w.min w.max w.var w.skew Dw.mean Dw.min Dw.max Dw.var Dw.skew]
% nvar=15;
% Z.mean=zeros(nmin*nhour*nday,nheight);
% Z.min =zeros(nmin*nhour*nday,nheight);
% Z.max =zeros(nmin*nhour*nday,nheight);
% Z.var =zeros(nmin*nhour*nday,nheight);
% Z.skew=zeros(nmin*nhour*nday,nheight);
% w.mean=zeros(nmin*nhour*nday,nheight);
% w.min =zeros(nmin*nhour*nday,nheight);
% w.max =zeros(nmin*nhour*nday,nheight);
% w.var =zeros(nmin*nhour*nday,nheight);
% w.skew=zeros(nmin*nhour*nday,nheight);
% Dw.mean=zeros(nmin*nhour*nday,nheight);
% Dw.min =zeros(nmin*nhour*nday,nheight);
% Dw.max =zeros(nmin*nhour*nday,nheight);
% Dw.var =zeros(nmin*nhour*nday,nheight);
% Dw.skew=zeros(nmin*nhour*nday,nheight);
% A=zeros(nheight*nvar,1);
% B=zeros(nheight,nvar);
% statfile=[way_proc_data_wband '1min_stat/VOCALS2008_1min_stat.dat'];
% fs=fopen(statfile,'r','ieee-le');
% frewind(fs);
ntime=32301; % number of minutes of data stored
% for iread=1:ntime
%     [t,count]=fread(fs,1,'int32');
%     if count>0                   % record found
%         time_aggreg(iread,1)=t;
%     else                         % catch and break for empty record
%         sprintf('iread=%d count=%d',iread, count)
%         break
%     end
% %     [kongflag(iread,1),count]=fread(fs,1,'int8'); % V6 2010-07-13
%     [angle,count]=fread(fs,4,'float32');
%     zen_angle.mean(iread)=angle(1);
%     zen_angle.std(iread) =angle(2);
%     zen_angle.cos(iread) =angle(3);
%     zen_angle.sin(iread) =angle(4);
%     %A=[Z.mean Z.min Z.max Z.var Z.skew w.mean w.min w.max w.var w.skew Dw.mean Dw.min Dw.max Dw.var Dw.skew]
%     [A,count]=fread(fs,nheight*nvar,'float64');
%     B=reshape(A,[nheight nvar]);
%     Z.mean(iread,:)=B(:,1); % 1st column
%     Z.min (iread,:)=B(:,2);
%     Z.max (iread,:)=B(:,3);
%     Z.var (iread,:)=B(:,4);
%     Z.skew(iread,:)=B(:,5);
%     w.mean(iread,:)=B(:,6);
%     w.min (iread,:)=B(:,7);
%     w.max (iread,:)=B(:,8);
%     w.var (iread,:)=B(:,9);
%     w.skew(iread,:)=B(:,10);
%     Dw.mean(iread,:)=B(:,11);
%     Dw.min (iread,:)=B(:,12);
%     Dw.max (iread,:)=B(:,13);
%     Dw.var (iread,:)=B(:,14);
%     Dw.skew(iread,:)=B(:,15);
% 
% end
% if iread~=ntime
%    sprintf('Warning: iread~=ntime')
% end
% 
% fclose(fs);
% truncate
% time_aggreg=time_aggreg(1:iread,:);
% kongflag=kongflag(1:iread,:); % V6 2010-07-13
% Z.mean=Z.mean(1:iread,:);
% Z.min =Z.min (1:iread,:);
% Z.max =Z.max (1:iread,:);
% Z.var =Z.var (1:iread,:);
% Z.skew=Z.skew(1:iread,:);
% w.mean=w.mean(1:iread,:);
% w.min =w.min (1:iread,:);
% w.max =w.max (1:iread,:);
% w.var =w.var (1:iread,:);
% w.skew=w.skew(1:iread,:);
% Dw.mean=Dw.mean(1:iread,:);
% Dw.min =Dw.min (1:iread,:);
% Dw.max =Dw.max (1:iread,:);
% Dw.var =Dw.var (1:iread,:);
% Dw.skew=Dw.skew(1:iread,:);

%cfad
% Z.cfad=NaN+zeros(ntime,nheight,nZbin,'uint8');
% w.cfad=NaN+zeros(ntime,nheight,nwbin,'uint8');
Dw.cfad=NaN+zeros(ntime,nheight,nDwbin,'uint8');
% AZ=zeros(nZbin*nheight,'uint8');
% Aw=zeros(nwbin*nheight,'uint8');
Ad=zeros(nDwbin*nheight,'uint8');
% BZ=zeros(nheight,nZbin,'uint8');
% Bw=zeros(nheight,nwbin,'uint8');
Bd=zeros(nheight,nDwbin,'uint8');
cfadfile=[way_proc_data_wband '1min_stat/onlywidth/VOCALS2008_cfad.dat'];
fc=fopen(cfadfile,'r','ieee-le');
frewind(fc);
for iread=1:ntime
%     [AZ,count]=fread(fc,nZbin*nheight,'uint8');
%     if count>0
%         BZ=reshape(AZ,[nZbin,nheight])';
%         Z.cfad(iread,:,:)=BZ; % (time, height, bin) = (time, 75, 38)
%     else
%         sprintf('%d',count)
%         iread=iread-1;
%         break
%     end
%     Aw=fread(fc,nwbin*nheight,'uint8');
%     Bw=reshape(Aw,[nwbin,nheight])';
%     w.cfad(iread,:,:)=Bw; % (time, height, bin) = (time, 75, 53)
    
    % Doppler width (fixed read 6/2011)
    Ad=fread(fc,nDwbin*nheight,'uint8');
    Bd=reshape(Ad,[nDwbin,nheight])';
    Dw.cfad(iread,:,:)=Bd; % (time, height, bin) = (time, 75, 73)
end
fclose(fc);

% time_aggreg in seconds since base_time
base_time=1225916508;
time_yday=datenum(0,0,0,0,0,base_time+time_aggreg)+datenum(1970,1,1,0,0,0)-datenum(2008,1,0,0,0,0); %yearday

% Z.bins=Zbins;
% Z.height=h;
% Z.time_yday=time_yday;
% w.bins=wbins;
% w.height=h;
% w.time_yday=time_yday;
Dw.bins=Dwbins;
Dw.height=h;
Dw.time_yday=time_yday;

% save as MATLAB variables
% save([way_proc_data_wband '1min_stat/Z_1min.mat'], 'Z','time_aggreg','time_yday');
% save([way_proc_data_wband '1min_stat/w_1min.mat'], 'w','time_aggreg','time_yday');
save([way_proc_data_wband '1min_stat/Dw_1min.mat'], 'Dw','time_aggreg','time_yday'); % Doppler width in v.6 added June 2011

% test Doppler width:
% imagesc(Dw.bins(2:end-1),Dw.height/1e3,squeeze(Dw.cfad(4501,:,2:end-1))); axis xy; axis([0 3 0 1.8])

% A netcdf file is made for cloud top height and cloud fraction 10-min statistics by proc_wband_cloudtop_10min.m