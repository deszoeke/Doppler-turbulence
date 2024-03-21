% cband_Z_histogram.m
%
% Make a histogram of the C-band column maximum reflectivity
% in 1 dBZ bins from 2-level PPI scans. The counts and relative probability in each bin
% have been computed by Casey Burleyson, NCSU. Sandra Yuter provides the range-dependent
% empirical noise/sensitivity, signals below which were excluded from the histogram.
% For lower reflectivity the radius of the area of sensitivity decreases. For C-band
% PPI probabilities I normalize the (pixel counts[Z])*(pixel area) by total area capable
% of sampling the reflectivity of each bin.
%
% Plot the probability computed from areal fraction of 
%
% Simon de Szoeke

read_parameters
domain_R=67; % km

B=load('~/Data/cruises/VOCALS_2008/RHB/radar/cband/C-Band_Column_Max_Reflectivity_Histogram.txt'); % Burleyson matrix
Zbin=B(:,1);
ccount=B(:,2);
frac=B(:,3);

% empirical noise level from Sandra
dBZ0=-33;
% Zmindetect=dBZ0+20*log10(r)

N=sum(ccount); % total number of pixels for entire cruise
Nppi=16410; % number of PPI scans included in histogram (Casey Burleyson)
apix=0.25*0.25; % pixel area, 0.25x0.25 km^2 (Casey Burleyson)
A=pi*domain_R^2; % total domain area

% compute the frequency as fractional area covered by the bin CMR
% divided by area in which this CMR is above noise
area_z=min(A,pi*10.^((Zbin-dBZ0)/10)); % maximum possible area (km^2) within which a signal equal to Zbin can be detected
% subtract out inner disk?
freq=(apix*ccount)./(area_z*Nppi); %  unitless area fraction
cfreq=1-cumsum(freq);
% note sum(freq)=0.50 of effectively-sampled area is counted in this distribution
%  N/Nppi*apix/A=0.07 is the fraction of the total domain area sampled

% %renormalize so that global integral of the probability distribution is 1
% freqn=freq./sum(freq);

cbz=[Zbin Zbin+1]';
cbf=[freq freq]';
cbz=cbz(:);
cbf=cbf(:);

% W-band CMR analysis
load([way_proc_data_wband '1min_stat_V5/Z_1min.mat']);
% define a noise level for the reflectivity based on Ken Moran's manual
radar_const=19.613; % dB_
dB_offset=-1.72; % dB
min_detectable_signal=-113.304+dB_offset; % dBmW
analog_noise=-120.74;  % peak, dBmW
digital_noise=-115.4;   % peak, dBmW
threshold_margin=3.5; % dB
noisefloor=20*log10(Z.height)+radar_const+digital_noise+threshold_margin;
% histogram of 1-minute vertical max dBZ (CMR) from W-band
[zmx,k]=max(Z.mean,[],2);
edges=-40:39;
wcount=histc(zmx,edges);

wbz=[edges' edges'+1]';
wbf=[wcount wcount]'/sum(wcount);
wbz=wbz(:);
wbf=wbf(:);

plot(Zbin(Zbin<4)+0.5, freqn(Zbin<4),'ro','markersize',2.5)
hold on
ii=find(cbz>=4);
plot(cbz(ii(2:end)),cbf(ii(2:end)),'r-')

ii=find(wbz>=-28);
plot(wbz(ii(2:end)),wbf(ii(2:end)),'b');
plot(edges(edges<-28)+0.5,wcount(edges<-28)/sum(wcount),'bo','markersize',2.5)
set(gca,'yscale','log','ylim',[1e-5 1e-1],'xlim',[-35 50],'fontsize',16)
ylabel('probability density p(Z<=CMR<Z+dZ) (unitless, dZ=1 dBZ)')
xlabel('reflectivity Z (dBZ)')
text(30,3e-4,'C-band','color','r','fontsize',14)
text(-20,3e-3,'W-band','color','b','fontsize',14)

print('-depsc',[way_proc_images_wband 'W+Cband_histogram_CMR.eps'])
print('-dpng',[way_proc_images_wband 'W+Cband_histogram_CMR.png'])