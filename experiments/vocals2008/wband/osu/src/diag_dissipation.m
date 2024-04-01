%diag_dissipation
% for debugging the dissipaton calculation step-by-step

factr=crosscomp/kolmogorov;

subplot(2,1,1,'align')
pcolor(-w'); shading flat; caxis([-1 1]);
hold on
plot(sum(isfinite(w)),1:120,'.-')
set(gca,'ylim',[40 55]);

subplot(2,1,2)
plot(-w(:,50:51));
axis tight; set(gca,'ylim',[-1 1])
hold on
% w visually correlated between top and next 1-2 levels

jj=isfinite(sum(w(:,50:51),2));
% sum(jj) % 4275, 4200 points is 20 minutes at 3.5 Hz
% sum(isfinite(w(:,50:51)))
% corr(w(jj,50),w(jj,51)); % corr at r=0.60, pretty good

%average top and 2d top when available
ser=nanmean(w(:,50:51),2);

% this is a procedure
% Identify the layer with the most coverage nonnoise returns
max(sum(isfinite(w)))
% identify the highest layer with 90% of the max coverage - cloud top is here or higher

% load 10-min ship-relative wind at cloud top interpolated from soundings
load('~/Data/cruises/VOCALS_2008/RHB/Scientific_analysis/programs/VOCALS2008_programs/wband/uv_cloudtop10.mat');
U=sqrt(utop10.^2+vtop10.^2); % m/s interpolated relative U from soundings and ship velocity
%yday_wband10=iday+idhr/24; % not used

%divide into 10-minute chunks
% filtering parameters
nwin=64;
F=(0:nwin/2)'/nwin*fs; % frequencies, Hz; 33: first is zero, last is Nyquist
dF=F(2)-F(1); % scalar
F53=F.^(5/3);
keep=2:8; % indices of spectra to keep
cut=10; % indices to regard as noise
dt=10*60; % seconds per window (10 min})
di=floor(fs*dt); % samples per window
for ind=1:round(nt/di) % index of the 10 minute interval in the hour
    i10=(1:di)+(ind-1)*di; % indices in the 10 minute interval
    ydayw10=iday+(idhr+(ind-1)/fs/3600)/24; % current time of this 10-minute window
    % Utop=nanmean(U(yday10>=ydayw10 & yday10<ydayw10+10/60/24));
    Utop=interp1(yday10,U,ydayw10,'nearest'); % interp ship-rel wind speed to the current time
    % what is the max coverage of nonnoise returns? --> nmax
    ii10=i10:min(i10+di-1,nt);
    if length(ii10)<di/2; continue; end % break from loop if too few obs.
    nw=sum(isfinite(w(ii10,:)));
    [nmax,knmax]=max(nw);
    % choose the highest level that has >0.5*nmax of this coverage --> ik
    ik=find(nw>0.5*nmax,1,'last');
    w10=w(i10:min(i10+di-1,nt),ik-1:ik+1);
    iwt=isfinite(w10);
    % take highest finite w within a level of this 50%ile cloud top
    wtop=w10(:,1); % w is positive-down
    wtop(iwt(:,2))=w10(iwt(:,2),2);
    wtop(iwt(:,3))=w10(iwt(:,3),3);
    % filter out extreme values
    wtop(abs(wtop-nanmean(wtop))>3*nanstd(wtop))=NaN;
    % interpolate missing values
    isn=isnan(wtop);
    wtop(isn)=interp1(time_yday_use(ii10(~isn)),wtop(~isn),time_yday_use(isn));
    [S,F]=pwelch(wtop(isfinite(wtop)),nwin,floor(nwin/2),nwin,fs); % uses hamming window 0.54-0.46cos
    %[Sm,Fm]=pwelch(wtop(isfinite(wtop)),hann(nwin),floor(nwin/2),nwin,fs); % hanning cos window (makes no difference)
    % low power in Nyquist frequency (end) could be from sampling out of phase with cycles
    Snoise=min(S(10:end-1)); % Fairall method, do not use Nyquist frequency
    vls=factr.*(2*pi/Utop)^(2/3).*mean(F53(keep).*(S(keep)-Snoise)); % dissipation ^2/3
    epsilon(ind,1)=vls.^1.5; % dissipation
    sigv2(ind,1)=dF*sum(S(2:end)-Snoise); % calculate variance from the resolved nonnoise spectrum
    %loglog(F,S,'.-','markersize',4)
end

% Farrar and Zappa suggest least squares fit for m and b to the linear spectrum S=m*x+b
% with x=k^-5/3, Snoise estimated by b. I worry least squares overfits to longer wavelengths
% due to the power law dependence.

% test plots
slo=mean(F53(keep).*(S(keep)-Snoise));
plot(1./F53,S,'.') % S vs. k-5/3
hold on
plot(1./F53,Snoise+slo./F53,'r') % Fairall model fit
set(gca,'yscale','log','xscale','log')
xlabel('k^{-5/3}')
ylabel('power spectral density')
% slope proportional to epsilon
% y-intercept is noise spectrum

plot(di*(1:nt/di)-di/2,3*sqrt(sigv2),'r-x',di*(1:nt/di)-di/2,1e3*epsilon,'r-o')
set(gca,'ylim',[-1 1.4])

clf
loglog(F,S,'m.-')
hold on
plot([1e-2 1e1],[1e1 1e-4],'r') % k^-5/3
axis([5e-2 2e0 1e-3 1e0])
set(gca,'fontsize',16)
xlabel('frequency (Hz)')
ylabel('power spectral density')