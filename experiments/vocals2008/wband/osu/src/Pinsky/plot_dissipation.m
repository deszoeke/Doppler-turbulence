% Calculate cloud top TKE dissipation from Pinsky-adjusted Doppler velocity spectrum
% Simon de Szoeke 2013.04.26

%% set environment and constants
%cd /Users/sdeszoek/Data/cruises/VOCALS_2008/RHB/Scientific_analysis/programs/VOCALS2008_programs/wband/Pinsky
warning('off','MATLAB:interp1:NaNinY')
run('../../read_parameters')
% way_Pinsky_retrieval='/Users/sdeszoek/Data/cruises/VOCALS_2008/RHB/Scientific_analysis/programs/VOCALS2008_programs/wband/Pinsky/retrieval/';
% way_Pinsky_retrieval='/Users/sdeszoek/Data/cruises/VOCALS_2008/RHB/Scientific_analysis/programs/wband/Pinsky/retrieval/';
way_Pinsky_retrieval='/Volumes/Pegasus2/cruises/VOCALS_2008/RHB/Scientific_analysis/programs/VOCALS2008_programs/wband/Pinsky';

yyyy=sprintf('%04d',year);
% stday=318; % why missing time adjusted motion till 318? % 281;
% enday=336;
stday=318;
enday=335;
fs=3.5; %Hz
crosscomp=3/4;
kolmogorov=0.54;
factr=crosscomp/kolmogorov;

%% structure function parameters
dxapprox=[1e3 5e2 1e2 5e1 1e1]; % m 
ndx=length(dxapprox);
% universal constant for 2nd order structure function
% longitudinal structure function constant C2 from Pope 6.2 (p.193) after Saddoughi and Veeravalli (1994)
C2ll=2.0;
factrz=1/C2ll;
factrx=3/4/C2ll;

% Kolmogorov (viscous) scale for atmospheric boundary layer is
% eta=(nu^3/epsilon); epsilon~3e-4 m^2/s^2/s, kinematic viscosity nu=2e-5 m^2/s
% --> eta= 2.3 mm
%     25 m = 1000 eta;  1000 m = 440 000 eta

%% w spectra window and FFT parameters used for dissipation
%nwin=64;
%keep=2:8; % indices of spectra to keep - not good if detrended
%cut=10; % index at which to start regarding as noise
% nwin=128;
% keep=4:16; % indices of spectra to keep, start at 4 for detrended & windowed spectra, ?32 eats into noise floor?
% cut=20; % index at which to start regarding as noise
nwin=256;
keep=8:32; % indices of spectra to keep, start at 4 for detrended & windowed spectra, ?32 eats into noise floor?
cut=40; % index at which to start regarding as noise


F=(1:nwin/2)'/nwin*fs; % frequencies, Hz; 64: first is nonzero, last is Nyquist
dF=F(2)-F(1); % scalar
F53=F.^(5/3);
dt=10*60; % seconds per window (10 min})
di=floor(fs*dt); % samples per window
S=NaN(size(F));

%% load 10-min ship-relative wind at cloud top interpolated from soundings
stem='';
load('/Volumes/Pegasus2/cruises/VOCALS_2008/RHB/Scientific_analysis/programs/VOCALS2008_programs/wband/uv_cloudtop10.mat');
load('/Volumes/Pegasus2/cruises/VOCALS_2008/RHB/Scientific_analysis/programs/VOCALS2008_programs/wband/uv10radar.mat');
% load('~/Data/cruises/VOCALS_2008/RHB/Scientific_analysis/programs/VOCALS2008_programs/wband/uv_cloudtop10.mat');
% load('~/Data/cruises/VOCALS_2008/RHB/Scientific_analysis/programs/VOCALS2008_programs/wband/uv10radar.mat');
% load('~/Data/cruises/VOCALS_2008/RHB/Scientific_analysis/programs/wband/uv_cloudtop10.mat');
% load('~/Data/cruises/VOCALS_2008/RHB/Scientific_analysis/programs/wband/uv10radar.mat');
% ship maneuvers removed
% horizontal ship-relative speed
U=sqrt(utoprel.^2+vtoprel.^2); % m/s interpolated relative U from soundings and ship velocity
Uradar=sqrt(uradarrel.^2+vradarrel.^2);

%% load data
load('cloudtopturb.mat')
% old version cloudtopturb_nfft128_2013-05-07.mat used nfft=128 points

% vertically-resolved data
load('radarturb.mat')
v2=load('radarturb2.mat')
DZ=(h(100)-h(99))*(1:5);
% old version radarturb_100m_2013-05-07.mat used only 100m structure function

%% recover noise for each spectrum
snoise=min(spec(:,cut:end-1),[],2);
snoisek=min(speck(:,cut:end-1),[],2);

%% QC for enough points in the spectrum
isenuf=count>1800;
isenufk=countk>1800;

% Already QCd Pinsky retrieval for missing motion data and motion
% correction failure, (85% good) and restored some motion values using
% the Crossbow accelerometer/gyro.
%  and, a posteriori, variance <0.4 m^2/s^2
isgdtmp=isenuf;
% nanstd(sigv2(isgdtmp))*3 ~ 0.4
isgd=isgdtmp & sigv2<0.4;
sp=spec(isgd,:);
ep=epsilon(isgd);
yd=ydayw10(isgd);

spk=speck(isgd,:);
epk=epsilonk(isgd,:);
% used for histogram

%%% plot
%% vertical structure of epsilon in cloud
clf
subplot(3,1,1)
hsol=area(ydayw10(isfinite(ydayw10))-304,max(0,-cos(2*pi*(ydayw10(isfinite(ydayw10))-5/24))),0,'edgecolor','w','facecolor','w');
set(gca,'color',0.8*[1 1 1])
hold on
% pcolor(ydayw10-304,h(1:100)/1e3,log10(epsx(:,:,3)'));shading flat
pcolor(ydayw10-304,h(1:100)/1e3,log10(epsilonk'));shading flat; caxis([-5 -2.5])
% pcolor(ydayw10-304,h(1:100)/1e3,log10(v2.epsilonk'));shading flat; caxis([-5 -2.5])
% pcolor(ydayw10-304,h(1:100)/1e3,(v2.epsilonk-epsilonk)');shading flat; caxis(1e-3*[-1 1])
hc=colorbar; set(hc,'fontsize',14)
set(hsol,'edgecolor','w','facecolor','w','visible','on')
set(get(hsol,'children'),'edgecolor','w','facecolor','w')
set(gca,'xminortick','on')
ylim([0 2])
set(gca,'fontsize',14,'tickdir','out','color',0.8*[1 1 1])
ylabel('height (km)')
xlabel('November 2008')
title('structure function dissipation log_{10}(\epsilon/[m^2 s^{-3}])')
orient landscape
box on
set(gca,'xlim',[15 32],'xtick',15:32)
set(gcf,'inverthardcopy','off','color','w')
% print -depsc epsilon_strfcn_timeheight.eps
% print -dpng epsilon_strfcn_timeheight.png
title('spectral dissipation [m^2 s^{-3}]')
% title('v2-v1 spectral dissipation [m^2 s^{-3}]')
print -depsc eps_spectr_v2_timeheight.eps
print -dpng eps_spectr_v2_timeheight.png

%% compare v2 epsilon with v1
clf
loglog(epsilonk(:),v2.epsilonk(:),'.','markersize',0.3)
axis square
axis([1e-6 1e-1 1e-6 1e-1])
hold on
plot([1e-6 1e-1],[1e-6 1e-1],'k-')
ylabel('v2'); xlabel('v1')
title('dissipation [m^2 s^{-3}]')
saveas(gcf,'compare_v1-v2_eps.png','png')

%% vertical top of dissipation meas.
isf=isfinite(epsx(:,:,3));
jmax=zeros(size(isf,1),1);
for j=1:size(isf,2)
    jmax(isf(:,j))=j;
end
hmax(jmax>0)=h(jmax(jmax>0));
hmax(jmax==0)=NaN;

% composite cloud-top relative diurnal epsilon
[epstmp,eps2tmp]=deal(zeros(144,100));
cnt=zeros(144,100);
tbin=(0:10:60*24)/60/24; % diurnal 10-min bins
isf=isfinite(ydayw10) & any(isfinite(epsx(:,:,3)),2);
fisf=find(isf);
iit=floor(interp1(tbin,1:length(tbin),mod(ydayw10(isf),1),'nearest'));
for it=1:length(iit)
    tmp=epsx(fisf(it),max(1,jmax(fisf(it))+(-99:0)),3);
    ii=isfinite(tmp);
    tmp(~ii)=0;
    epstmp(iit(it),:)=epstmp(iit(it),:)+tmp;
    eps2tmp(iit(it),:)=eps2tmp(iit(it),:)+tmp.*tmp;
    cnt(iit(it),:)=cnt(iit(it),:)+ii;
end
epscmp=epstmp./cnt;
epscmp(cnt<1)=NaN;
epsstd=sqrt(eps2tmp./cnt-epscmp.*epscmp);
epsstd(cnt<4)=NaN;

clf
subplot(2,1,1)
hsol=area(tbin*24,-2+0.4*max(0,-cos(2*pi*(tbin-5/24))),-2,'edgecolor','w','facecolor','w');
hold on
pcolor(tbin(1:end)*24,(h(19:100)-h(100))/1e3,log10(epscmp([1:end 1],19:end))'); shading flat
colorbar
axis([0 24 -2 0])
caxis([-5 -2.5])
set(gca,'color',0.8*[1 1 1],'fontsize',14,'xtick',0:4:24,'tickdir','out')
set(hsol,'edgecolor','none')
ylabel('height below cloud top (km)')
title('log10 dissipation [m^2 s^{-3}]')

subplot(2,1,2)
hsol=area(tbin*24,-2+0.4*max(0,-cos(2*pi*(tbin-5/24))),-2,'edgecolor','w','facecolor','w');
hold on
pcolor(tbin(1:end)*24,(h(19:100)-h(100))/1e3,(epsstd([1:end 1],19:end)./epscmp([1:end 1],19:end))'); shading flat
colorbar
axis([0 24 -2 0])
caxis([0 2])
set(gca,'color',0.8*[1 1 1],'fontsize',14,'xtick',0:4:24,'tickdir','out')
set(hsol,'edgecolor','none')
title('<dissipation''^2>^{1/2}/dissipation')
xlabel('hour UTC')

set(gcf,'inverthardcopy','off','color','w')
orient tall
print -depsc dissipation_diurnal_cmp.eps

% Standard deviation of epsilon is on the order of
% epsilon itself, indicating to me turbulence is not terribly episodic.

% note inversion properties are saved in
% /Users/sdeszoek/Data/cruises/VOCALS_2008/RHB/balloon/Processed/sonde_inversion.mat