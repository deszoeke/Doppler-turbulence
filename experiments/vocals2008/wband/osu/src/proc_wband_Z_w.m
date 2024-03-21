% proc_wband_Z_w.m
% 2009-04-08 :: VOCALS 2008 :: Simon de Szoeke
%
% Simon's first procedural attempt to process the W-band radar data for one
% file/hour:
% 1. read radar data
% 2. read motion correction file
% 3. check filtered velocities vs. ship motion
% 4. Subtract ship motion from w.
% 5. Plot reflectivity and W

%cd /Users/sdeszoek/Data/cruises/VOCALS_2008/RHB/Scientific_analysis/programs/VOCALS2008_programs/wband
%addpath('~/Data/cruises/VOCALS_2008/RHB/Scientific_analysis/programs/VOCALS2008_programs/');
read_parameters;

year=2008;
yday=326;
hour=4;
% d326 h04 is a good time to study the noise --very little signal
% d330 h06 the motion correction was not working; roll angle parked at -8 degrees
yyyy=sprintf('%04d',year);
ddd= sprintf('%03d',yday);
hh=  sprintf('%02d',hour);

% Read radar data
momentfile=dir([way_raw_data_wband yyyy ddd hh '*MMCRMom.nc']);
if length(momentfile)>1
    disp(['Multiple files found: ' momentfile(:).name])
end
filename=[way_raw_data_wband momentfile.name];
mm=momentfile.name(10:11);
%nc_dump(filename)

% base_time in seconds since 1970-1-1 00:00:00
base_time=nc_varget(filename,'base_time');
% base_time_mld in matlab datenumber
base_time_mld=double(base_time)/86400 + datenum(1970,1,1,0,0,0);
base_time_yday=base_time_mld-datenum(year,1,0,0,0,0);
% time_offset in seconds
time_offset=nc_varget(filename,'time_offset');

range=nc_varget_lapxm(filename,'Heights',[0 0],[1 -1]);
% define a noise level for the reflectivity based on Ken Moran's manual
dB_offset=-1.72; % dB, Ken Moran September 2009 recalibration
radar_const=19.613; % dB_
min_detectable_signal=-113.304+dB_offset; % dBmW
analog_noise_signal=-120.74; % dBmW peak diagnosed by Simon
digital_noise_signal=-115.4;
noise_margin=3.5; % dB
% min_detectable_signal is good basis for a threshold, not to offset
% meteorological signal for the noise.
noisefloor=20*log10(range)+radar_const+digital_noise_signal+noise_margin; % dBZ
noiseavg=20*log10(range)+radar_const+analog_noise_signal;

%threshold=nc_varget_lapxm(filename,'MinimumDetectableReflectivity',[0 0],[1 -1]);
RCpower=nc_varget_lapxm(filename,'RangeCorrectedPower') +dB_offset;
Z=nc_varget_lapxm(filename,'Reflectivity')              +dB_offset;
w=nc_varget_lapxm(filename,'MeanDopplerVelocity');
snr=nc_varget_lapxm(filename,'SignalToNoiseRatio');
noise=nc_varget_lapxm(filename,'NoiseLevel')            +dB_offset;

Zrec=10.^(Z/10);
Zmet=Zrec - repmat(10.^(noiseavg/10),[size(Z,1) 1]);
dBZmet=10*log10(Zmet);
% This does not account for the finite width of the noise, which will broaden the
% reflectivity signals. I should still threshold the data below the noisefloor

% Read motion-correction files with Alan Brewer's time adjustment
kongfile=dir([way_raw_data_wband 'motion_adjT/' yyyy ddd hh '*Kongsberg_adjT.txt']);
if isempty(kongfile)
    disp('No Kongsberg file in hour. Files in previous hour:')
    ls([way_raw_data_wband 'motion_adjT/' yyyy ddd sprintf('%02i',hour-1) '*Kongsberg_adjT.txt'])
elseif length(kongfile)>1
    disp(['Multiple files found: ' kongfile(:).name])
end
%kongfilename=[way_raw_data_wband 'motion/20083300600Kongsberg.txt']; % not time-shifted file
kongfilename=[way_raw_data_wband 'motion_adjT/' kongfile(1).name];
[kongtime,kongw,kongpitch,kongroll]=read_kongsberg(kongfilename);

% synchronize Kongsberg and radar time (yearday)
kongw4radar=interp1(kongtime,kongw,base_time_yday+time_offset/86400);
pitch=interp1(kongtime,kongpitch,base_time_yday+time_offset/86400);
roll=interp1(kongtime,kongroll,base_time_yday+time_offset/86400);
% motion-compensated drop vertical velocity +up
wdrop=-w-repmat(kongw4radar,[1,120]);

% pitch/roll corrected vertical coordinate
quad=sin(pitch/180*pi).^2+sin(roll/180*pi).^2; % equals sin(theta)^2
sintheta=sqrt(quad);
theta=asin(sintheta); % zenith angle (radians)
costheta=sqrt(1-quad);
height=costheta * range;

% correlate Kongsberg motion with top of cloud Doppler w
icltop=range>1500 & range<1600;
%plot(kongtime,kongw,base_time_yday+time_offset/86400,kongw4radar,'r.')
C=nancov(kongw4radar,-mean(w(:,icltop),2));
corw=C(2,1)/sqrt(C(1,1)*C(2,2));

clf
plot(kongw4radar,-mean(w(:,icltop),2),'b.','markersize',2)
axis equal
hold on
plot([-2 2],[-2 2],'k')
xlabel('Kongsberg vertical velocity (m s^{-1})')
ylabel('W-band Doppler vertical velocity (m s^{-1})')
text(0.5,-2,sprintf('r=%4.2g',corw));

% PLOT
mask=Z>repmat(noisefloor,[length(time_offset) 1]);
wdropmask=wdrop;
wdropmask(~mask)=NaN;
b2rcolormap(16);
subplot(2,1,1)
imagesc(time_offset/60,range(1,4:end)',Z(:,4:end)');
set(gca,'ydir','normal')
hc=colorbar; set(get(hc,'ylabel'),'string','dB')
%xlabel([datestr(base_time_mld,'yyyy-mmm-dd HH:MM') ' time (minutes)'])
ylabel('height (m)')
title({momentfile.name 'reflectivity (dBZ)'})
subplot(2,1,2)
pcolor(time_offset/60,range(1,4:end),double(wdropmask(:,4:end))')
shading flat
caxis([-4 4])
set(gca,'ydir','normal','color',0.7+[0 0 0])
hc=colorbar; set(get(hc,'ylabel'),'string','m s^{-1}')
xlabel([datestr(base_time_mld,'yyyy-mmm-dd HH:MM') ' time (minutes)'])
ylabel('height (m)')
title('motion-compensated (+up) vertical velocity (m s^{-1})')
%print -depsc2 SNR_W.eps

if false % Plot Reflectivity (signal) and Noise
b2rcolormap(16);
subplot(5,1,1)
imagesc(time_offset/60,height(1,4:end),Z(:,4:end)')
set(gca,'ydir','normal')
hc=colorbar; set(get(hc,'ylabel'),'string','dBZ')
%xlabel(datestr(base_time_mld,'yyyy-mmm-dd HH:MM') ' time (minutes)'])
ylabel('height (m)')
title({momentfile.name ' Reflectivity (dBZ)'})
subplot(5,1,2)
pcolor(time_offset/60,height(1,4:end),noise(:,4:end)')
shading flat
caxis([60 100])
set(gca,'ydir','normal','color',[0 0 0])
hc=colorbar; set(get(hc,'ylabel'),'string','dB')
%xlabel([datestr(base_time_mld,'yyyy-mmm-dd HH:MM') ' time (minutes)'])
ylabel('height (m)')
title('Noise (dB)')
subplot(5,1,3)
pcolor(time_offset/60,height(1,4:end),(Z(:,4:end)-noise(:,4:end))')
shading flat
set(gca,'ydir','normal','color',[0 0 0])
hc=colorbar; set(get(hc,'ylabel'),'string','dB')
%xlabel([datestr(base_time_mld,'yyyy-mmm-dd HH:MM') ' time (minutes)'])
ylabel('height (m)')
title('Signal(dB)-Noise(dB)')
subplot(5,1,4)
pcolor(time_offset/60,height(1,4:end),snr(:,4:end)')
shading flat
set(gca,'ydir','normal','color',[0 0 0])
hc=colorbar; set(get(hc,'ylabel'),'string','dB')
%xlabel([datestr(base_time_mld,'yyyy-mmm-dd HH:MM') ' time (minutes)'])
ylabel('height (m)')
title('SNR (dB)')
subplot(5,1,5)
pcolor(time_offset/60,height(1,4:end),(Z(:,4:end)-noise(:,4:end)-snr(:,4:end))')
shading flat
set(gca,'ydir','normal','color',[0 0 0])
hc=colorbar; set(get(hc,'ylabel'),'string','dB')
xlabel([datestr(base_time_mld,'yyyy-mmm-dd HH:MM') ' time (minutes)'])
ylabel('height (m)')
title('Signal-Noise-SNR (dB)')
end

% CFAD of Z
bin=-65:.5:-20; % lower edges of bins
count=histc(Z,bin);

% plot CFAD
[cc,hc]=contourf(bin+(bin(2)-bin(1))/2,range,max(0,log10(count))',40);
set(hc,'edgecolor','none')
%pcolor(bin,height,count') % shows discrete ranges
%shading flat
%caxis([0 max(count(:))/4])
% axis([-65 0 0 3000])
% colormap([ones(2,3); b2rcolormap(28)]); % white for printing
colormap(b2rcolormap(40));
hb=colorbar('southoutside');
set(get(hb,'xlabel'),'string','count (decade)')
title({[yyyy ' ' ddd ' ' hh ':00'] 'Reflectivity CFAD'})
xlabel('dBZ')
ylabel('range (m)')
hold on
plot(noisefloor,height,'k--')
%print('-dpng',[way_proc_images_wband 'CFAD_Z20083260400_noise.png'])


% CFAD of range-UNcorrected Z
rucZ=Z-repmat(20*log10(range),[length(time_offset) 1]);
binrucz=-105:0.2:-85;
countrucz=histc(rucZ,binrucz);
[ch,hc]=contourf(binrucz+(binrucz(2)-binrucz(1))/2,range,max(0,log10(countrucz))',40);
set(hc,'edgecolor','none')

% CFAD of power
power=nc_varget_lapxm(filename,'Power');
binp=55:135;
countp=histc(power,binp);
binpnoise=55:0.01:65;
countp_norange=histc(power(:),binpnoise);
figure;dock;
[cc,hc]=contourf(binp+(binp(2)-binp(1))/2,range,max(0,log(countp))',40);
set(hc,'edgecolor','none')
%caxis([0 max(count)/2])
axis([55 75 0 3000])
colormap(b2rcolormap(30));
hb=colorbar('southoutside');
set(get(hb,'xlabel'),'string','log(count)')
title({[yyyy ' ' ddd ' ' hh ':00'] 'Power CFAD'})
xlabel('dB?')
ylabel('range (m)')
hold on
plot(binpnoise,countp_norange,'m.-');
orient tall
% peak of power distribution is 58.76 dBmW
% min. detectable signal dBmW = -113.304

% CFAD of range-UNcorrected range corrected power
binruc=-125:-40;
rucpower=RCpower-repmat(20*log10(range),[length(time_offset) 1]);
countruc=histc(rucpower,binruc);
binrucnoise=-122:.02:-100;
countnoise=histc(rucpower(height>100 & height<=3000),binrucnoise); % store/save as uint16
[cc,hc]=contourf(binruc+(binruc(2)-binruc(1))/2,range,max(0,log10(countruc))',40);
set(hc,'edgecolor','none')
hold on
plot(binrucnoise,log10(countnoise)*5e2,'m.-');
axis([-122 -100 0 3000])
plot(min_detectable_signal+[0 0],[0 3000],'k--')
ylabel('range(m)')
title({['Histogram ' yyyy ' ' ddd ' ' hh],'W-band calibrated range-UNcorrected power (dBmW)'})
hcb=colorbar('southoutside');
set(get(hcb,'xlabel'),'string','log_{10} frequency')
text(-112,2600,'min. detectable power','horizontalalignment','right')
ht=text(-116.5,700,{'power histogram', '500 m = 1 decade'},'color','m');
%print('-dpng',[way_proc_images_wband 'CFAD_P20083260400_noise.png'])

% CFAD of range-corrected power
binrcp=-84:1:4;
countrcp=histc(RCpower,binrcp);
binrcpnoise=55:0.01:65;
countp_norange=histc(power(:),binpnoise);
figure;dock;
%[cc,hc]=contourf(binrcp+(binrcp(2)-binrcp(1))/2,range,max(0,log(countrcp))',40);
[cc,hc]=contourf(bintest+(bintest(2)-bintest(1))/2,range,max(0,log(test))',40);
set(hc,'edgecolor','none')
%caxis([0 max(count)/2])
axis([-83 4 0 3000])
colormap(b2rcolormap(40));
hb=colorbar('southoutside');
set(get(hb,'xlabel'),'string','log(count)')
title({[yyyy ' ' ddd ' ' hh ':00'] 'Power CFAD'})
xlabel('dBmW')
ylabel('range (m)')
hold on
plot(binpnoise,countp_norange,'m.-');
orient tall
