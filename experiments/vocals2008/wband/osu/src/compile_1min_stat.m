% compile_wband_1min_stat.m
% VOCALS 2008 :: 2009-07-07 :: Simon de Szoeke
%
% Compile 5 vertically-resolved statistics of Z and w:
% mean, min, max, variance, and skewness
%
% Compile vertically-resolved histograms of reflectivity and w.
%
% Input serial unformatted binary files from proc_wband_1min_stat.m.
% Compile cloud top into a netcdf file with proc_wband_cloudtop_10min.m.
%
% Dimensions:
% time( )                               serial time (s), record dimension
% height(75)    134:1982                center of range gates (m)
% Zbins(38)     [-Inf -40:2:30 Inf]     reflectivity histogram bin edges (dBZ)
% wbins(53)     [-Inf -5:0.2:5 Inf]     velocity histogram bins edges (m/s)

read_parameters;

% concatenate 1-min data files if not done already
if ~exist([way_proc_data_wband '1min_stat/VOCALS2008_1min_stat.dat'],'file')
    filecat([way_proc_data_wband '1min_stat/2008*_1min_stat.dat'],...
            [way_proc_data_wband '1min_stat/VOCALS2008_1min_stat.dat'])
end
if ~exist([way_proc_data_wband '1min_stat/VOCALS2008_cfad.dat'],'file')
    filecat([way_proc_data_wband '1min_stat/2008*_cfad.dat'],...
            [way_proc_data_wband '1min_stat/VOCALS2008_cfad.dat'])
end

% Define dimensions
% dBZ and w bin edges
Zbins=[-inf -40:2:30 inf];
wbins=[-inf -5:0.2:5 inf];
% vertical coordinate (center of ranges, m)
h=load([way_proc_data_wband '1min_stat/height.txt']);
nheight=length(h); % 75
nday=25;
nmin=60;
nhour=24;
nZbin=length(Zbins);
nwbin=length(wbins);

%mean(time[ ],height[75])
time_aggreg=zeros(nmin*nhour*nday,1);
zen_angle.mean=NaN+zeros(nmin*nhour*nday,1);
zen_angle.std=zen_angle.mean;
zen_angle.cos=zen_angle.mean;
zen_angle.sin=zen_angle.mean;
Z.mean=zeros(nmin*nhour*nday,nheight);
Z.min =zeros(nmin*nhour*nday,nheight);
Z.max =zeros(nmin*nhour*nday,nheight);
Z.var =zeros(nmin*nhour*nday,nheight);
Z.skew=zeros(nmin*nhour*nday,nheight);
w.mean=zeros(nmin*nhour*nday,nheight);
w.min =zeros(nmin*nhour*nday,nheight);
w.max =zeros(nmin*nhour*nday,nheight);
w.var =zeros(nmin*nhour*nday,nheight);
w.skew=zeros(nmin*nhour*nday,nheight);
A=zeros(nheight*10,1);
B=zeros(nheight,10);
statfile=[way_proc_data_wband '1min_stat/VOCALS2008_1min_stat.dat'];
fs=fopen(statfile,'r','ieee-le');
frewind(fs);
iread=0;
while ~feof(fs)
    iread=iread+1;
    [t,count]=fread(fs,1,'int32');
    if count>0 % catch and break for empty record
        time_aggreg(iread,1)=t;
    else
        sprintf('%d',count)
        iread=iread-1;
        break
    end
    [angle,count]=fread(fs,4,'float32');
    zen_angle.mean(iread)=angle(1);
    zen_angle.std(iread) =angle(2);
    zen_angle.cos(iread) =angle(3);
    zen_angle.sin(iread) =angle(4);
    [A,count]=fread(fs,nheight*10,'float64');
    B=reshape(A,[nheight 10]);
    Z.mean(iread,:)=B(:,1); % 1st column
    Z.min (iread,:)=B(:,2);
    Z.max (iread,:)=B(:,3);
    Z.var (iread,:)=B(:,4);
    Z.skew(iread,:)=B(:,5);
    w.mean(iread,:)=B(:,6);
    w.min (iread,:)=B(:,7);
    w.max (iread,:)=B(:,8);
    w.var (iread,:)=B(:,9);
    w.skew(iread,:)=B(:,10);
end
fclose(fs);
% truncate
time_aggreg=time_aggreg(1:iread,:);
Z.mean=Z.mean(1:iread,:);
Z.min =Z.min (1:iread,:);
Z.max =Z.max (1:iread,:);
Z.var =Z.var (1:iread,:);
Z.skew=Z.skew(1:iread,:);
w.mean=w.mean(1:iread,:);
w.min =w.min (1:iread,:);
w.max =w.max (1:iread,:);
w.var =w.var (1:iread,:);
w.skew=w.skew(1:iread,:);
ntime=iread;

%cfad: [38x75x60 double]
Z.cfad=NaN+zeros(ntime,nheight,nZbin,'uint8');
w.cfad=NaN+zeros(ntime,nheight,nwbin,'uint8');
AZ=zeros(nZbin*nheight,'uint8');
Aw=zeros(nwbin*nheight,'uint8');
BZ=zeros(nheight,nZbin,'uint8');
Bw=zeros(nheight,nwbin,'uint8');
cfadfile=[way_proc_data_wband '1min_stat/VOCALS2008_cfad.dat'];
fc=fopen(cfadfile,'r','ieee-le');
frewind(fc);
for iread=1:ntime
    [AZ,count]=fread(fc,nZbin*nheight,'uint8');
    if count>0
        BZ=reshape(AZ,[nZbin,nheight])';
        Z.cfad(iread,:,:)=BZ; % (time, height, bin) = (time, 75, 38)
    else
        sprintf('%d',count)
        iread=iread-1;
        break
    end
    Aw=fread(fc,nwbin*nheight,'uint8');
    Bw=reshape(Aw,[nwbin,nheight])';
    w.cfad(iread,:,:)=Bw; % (time, height, bin) = (time, 75, 53)
end
fclose(fc);

% time_aggreg in seconds since base_time
base_time=1225916508;
time_yday=datenum(0,0,0,0,0,base_time+time_aggreg)+datenum(1970,1,1,0,0,0)-datenum(2008,1,0,0,0,0); %yearday

Z.bins=Zbins;
Z.height=h;
Z.time_yday=time_yday;
w.bins=wbins;
w.height=h;
w.time_yday=time_yday;
% save as MATLAB variables
save([way_proc_data_wband '1min_stat/Z_1min.mat'], 'Z','time_aggreg','time_yday');
save([way_proc_data_wband '1min_stat/w_1min.mat'], 'w','time_aggreg','time_yday');

% A netcdf file is made for cloud top height and cloud fraction 10-min statistics by proc_wband_cloudtop_10min.m