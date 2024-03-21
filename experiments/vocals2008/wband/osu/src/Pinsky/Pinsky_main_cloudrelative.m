% Load hourly radar data, do corrections, feed to Pinsky_retrieval

%% set environment and constants
%cd /Users/sdeszoek/Data/cruises/DYNAMO_2011/Revelle/Scientific_analysis/programs/VOCALS2008_programs/wband/Pinsky
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
% year=2011;
yyyy=sprintf('%04d',year);
stday=318;
enday=336;
fs=3.5; %Hz

%% parameters of Pinsky retrieval
%Zbin=-40:0.25:30;
Zbin=-50:0.5:30; % adequate and more robust than 0.25 dBZ bins
thickbin=0:100:1400; % SPdeS added cloud thickness bins
rho=0;
%rho=[0:64 repmat(64,1,55)]*.2/64; % approx correlation from Pinsky et al. Fig.5

% weights for a time-localized lookup table for Vg(hrel,Z)
% also used to filter 1-min cloud top height
iw=-12:12; % offsets, 0 is centered
weight=exp(-(iw/3).^2);

%% load 1 min cloud fraction and top height data
A=load([way_proc_data_wband 'cloudheight/CloudHeight_1min_2008310-336.txt']);
A(A==-999)=NaN; % NaN out missing value of -999
time=A(:,1);
cloud_fraction=A(:,2);
cloud_top_mean=A(:,3);
% cloud_top_std=A(:,4);
% cloud_top_val15=A(:,5);
% cloud_top_median=A(:,6);
% cloud_top_val85=A(:,7);
cloud_top_valmean=A(:,8); %mean of tops within 15-85 percentile
num_returns=A(:,9);
num_cloud_returns=A(:,10);
% num_clear_returns=A(:,11);
% num_nonmet_returns=A(:,12);
flag=A(:,13);
% flag values
% 0 no errors detected, cloud top computed. 
% 1 cloud fraction=0, cloud top set to NaN.
% 2 cloud fraction>=0 but found no 3 cloud returns with consistent height,
%   possible glitch, provisional cloud top computed.
% NaN out percentile-based cloud tops=0 using flag
lf=logical(flag);
% cloud_top_val15(lf)=NaN;
% cloud_top_val85(lf)=NaN;
cloud_top_valmean(lf)=NaN;
cloud_top_median(lf)=NaN;
cloudtop=spatialfilter(cloud_top_mean,weight,1,100,100); %filtered and filled
starter=1;
base_time(starter)=1225916508; % s since 1970-1-1 00:00
time_yday_top=datenum(0,0,0,0,0,1225916508+time)+datenum(1970,1,1,0,0,0)-datenum(2008,1,0,0,0,0);

%% begin loading Wband moment data independent height variable
momentfile=dir([way_raw_data_wband yyyy '*MMCRMom.nc' ]);
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
foot=0.01; % minimum Ustd (m/s) used to get S1,S2, reduces a rare numerical glitch
% hthick=1; % kludge, must be defined but not used

%% 1st huge loop over all radar files by day, hour, making+saving hourly Vg lookup tables
%if ~exist('./retrieval/Pinsky_lookup2008_335_23_cloudrelative_thickness.mat','file')
if ~exist('./retrieval/Pinsky_lookup2008_335_23_cloudrelative_p5dBZ.mat','file')
    Pinsky_loopread_cloudrelative_p5dBZ
end
hrel=h(1:101)-h(91);

%% 2nd huge loop 
Pinsky_loopread_retrieval_cloudrelative

return
