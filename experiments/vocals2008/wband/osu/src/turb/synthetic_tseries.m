% Create a synthetic 1D velocity time series for the inertial subrange.

% sampling parameters
N=2100;        % number of points in time series
fsampl=3.5;    % sampling frequency
fnyq=fsampl/2;
f=fsampl*(0:N-1)'/N;
Tsampl=1/fsampl;
Ttot=N*Tsampl;
t=(0:N-1)'*Tsampl;

omegasampl=2*pi*fsampl;
omega=omegasampl*(0:N-1)'/N;    % sample frequencies (rad/second)
omeganyq=omegasampl/2;

% convert to spatial coordinate and wavenumbers
U=10; % m/s mean velocity in the x direction
DXsampl=Tsampl*U;
Xtot=Ttot*U;  % Xtot=Lsampl  length of series (m)
x=(0:N-1)'*DXsampl; % 0-6 km

k=omega/U;    % sample wavenumbers (rad/meter)
ksampl=omegasampl/U;
knyq=omeganyq/U;

% create synthetic ideal series of length M>N from which x and k subsample
M=2^14;
DXideal=DXsampl/4;         % ideal mesh is finer
xideal=DXideal*(0:M-1)';
kidealup=2*pi*(0:M-1)'/(M*DXideal);
kideal=[kidealup(1:M/2+1); -kidealup(M/2:-1:2)]; % 2nd half negative wavenumbers
Dkideal=2*pi/(M*DXideal);
kidealsampl=2*pi/DXideal;
% highest ideal wavenumer is max(kideal)=2pi/(0.7m), wavelength=0.7 m, inertial
% highest resolved wavenumber is max(k)=2*pi/2.9m

% spectral energy density function
epsilon=.01;     % dissipation m2/s3
% spectral energy density const. for velocity transverse to sampling, e.g E_zz(k_x) (eqn 6.243 Pope)
C=0.65;
E=C*epsilon^(2/3)*abs(kideal).^(-5/3); % double counting for negative k
E(1)=0; % zero out first wavenumber, k=0
varE=sum(E)*Dkideal;

% complex spectral coefficients u(kideal) of velocity wideal(xideal)
phase=2*pi*rand(M,1);  % random phase shifts
phase(M/2+2:end)=-phase(M/2:-1:2); % symmetrize phase shifts to get real w

% 1/2 fixes double counting conjugates at negative wavenumbers
edens=E/2;
edens(1)=0;
edens(M/2+1)=E(M/2+1);
% integrate spectral density to spectral coefficient
u=sqrt(edens*Dkideal).*exp(1i*phase);
% velocity spatial series wideal(xideal)
wideal=real(fft(u));

% ideal spectrum at full resolution
[sideal,ktest]=pwelch(wideal, 1024, 512   , 1024,kidealsampl);
% does not have reduced power law "foot" at high wavenumber;
% presume that the foot is a result of resampling.

% add a spike of frequency kspike
kspike=0.1; % 1/m
A=1e3*E(find(kideal>kspike,1,'first'));
spikeideal=A*sin(kspike*xideal);
spike=A*sin(kspike*x);

% subsample w in pulses on x mesh
%w=interp1(xideal,wideal,x); % interpolate
Tdwell=0.123; % radar dwell time (s)
Xdwell=Tdwell*U;
w=dwellavg(xideal,wideal,x-Xdwell/2,Xdwell); % centered average dwell-length intervals
% plot(xideal,wideal)
% hold on
% plot(x,w,'r')

[P,  kwelch]=pwelch(w, 512, 256   , 512,ksampl);  % PSD (normalization?)
%                   window overlap  nfft
[P_s,kwelch_s]=specsmoo(P,ksampl);
L=length(P);
s=fft(w);
S=(s(1:N/2+1).*conj(s(1:N/2+1)))/N; % power spectrum correctly normalized
[S_s,k_s]=specsmoo(S,ksampl);

[Ppsd,kpsd]=psd2(w,512,ksampl,512,256);
[Ppsd_s,kpsd_s]=specsmoo(Ppsd,ksampl);
% specsmoo allows us to use an appropriate window for med-low wavenumbers
% but averages together adjacent spectral estimates from higher wavenumbers
% to reduce noise. The wavenumbers output are midpoints of the wavenumber
% averaging intervals. The width of the wavenumber intervals from specsmoo
% grows roughly exponentially.

% use Fairall wband_velsp4.m method
[S_w,K]=psd2(detrend(w),length(w),ksampl);%power spectrum
[S_ws,Ks]=specsmoo(S_w,ksampl);%Smoothed spectrum

% window design:
% N is 2100 for a 10-minute interval
nwindow=420; % 2 minute
noverlap=0.5*nwindow;
nws=9;       % 9 complete windows in the 10-minute interval
nfft=512;
% use spectrum objects to compute the spectrum
hwelch=spectrum.welch('Hann',nwindow,100*noverlap/nwindow);
hpsd=psd(hwelch,w,'NFFT',nfft,'Fs',ksampl);
[sdata_s,sfreq_s]=specsmoo(hpsd.data,ksampl);

%plot spectrum
clf
loglog(k_s,S_s,'r') % only to Nyquist wavenumber
hold on
loglog(kwelch_s,P_s,'b')
loglog(kpsd_s,Ppsd_s,'g')
loglog(sfreq_s,sdata_s,'m')
loglog(kideal,E,'k--')
legend('FFT','Welch','PSD2','spectrum','ideal')
title('smoothed spectra of a synthetic \itk^{-5/3}\rm series')
xlabel('wavenumber \itk\rm (m^{-1})')
ylabel('spectral energy density \itE\rm (m^2s^{-2})')
text(3e-1,4e-3,'E=C\epsilon^{-2/3}k^{-5/3}')
%print('-dpng','synthetic_k53spectra.png')

% compute epsilon from each wavenumber in the spectrum
eps_ideal=(E/C).^(3/2).*abs(kideal).^(5/2);
eps_welch=(P/C).^(3/2).*abs(kwelch).^(5/2);
eps_spect=(hpsd.data/C).^(3/2).*abs(hpsd.frequencies).^(5/2);
eps_psd=(Ppsd/C).^(3/2).*abs(kpsd).^(5/2);
eps_fft=(S/C).^(3/2).*abs(k(1:N/2+1)).^(5/2);

% figure; dock
clf
subplot(2,1,1)
plot(k(1:N/2+1),eps_fft,'r')
hold on
plot(kwelch,eps_welch,'b')
plot(kpsd,eps_psd,'g')
plot(hpsd.frequencies,eps_spect,'m')
plot(kideal(kideal>0),eps_ideal(kideal>0),'k--')
set(gca,'xlim',[0 1.2])
ylabel('dissipation \epsilon (m^2 s^{-3})')
xlabel('k')
subplot(2,1,2)
plot(k(1:N/2+1),cumsum(eps_fft)*k(2)./k(1:N/2+1),'r')
hold on
plot(kwelch,cumsum(eps_welch)*kwelch(2)./kwelch,'b')
plot(kpsd,cumsum(eps_psd)*kpsd(2)./kpsd,'g')
plot(hpsd.frequencies,cumsum(eps_spect)*hpsd.frequencies(2)./hpsd.frequencies,'m')
plot(kideal(kideal>0),cumsum(eps_ideal(kideal>0))*kideal(2)./kideal(kideal>0),'k--')
set(gca,'xlim',[0 1.2])
legend('FFT','Welch','PSD2','spectrum','ideal',0)
xlabel('wavenumber k (m^{-1})')
ylabel('k^{-1}\int_0^k\epsilon dk''')

% Integrate the energy in each spectrum in the range k(2:L/2+1) to find the normalization ratios
% 2*pi/k(2) = 6 km          (energy containing scale)
% 2*pi/k(L/2+1) = 23.4 m    (inertial)       -- L chosen by pwelch
% 2*pi/k(N/2+1) = 5.7 m
kii=kideal>=kwelch(2) & kideal<=kwelch((L+1)/2);
ki=k>=kwelch(2) & k<=kwelch((L+1)/2);
% mostly influenced by low wavenumbers
Sint=sum(S(2:(L+1)/2));
Pint=sum(P(2:(L+1)/2));
Eint=sum(E(kii));

% Windowing filters the variance let through the PSD, so compute the
% variance let through the windows
startwindow=1+(0:0.5*nwindow:floor(N-nwindow));
win1=zeros(N,1);
win1(1:nwindow)=window('Hann',nwindow);
win1(nwindow+1:N)=0;
win=zeros(N,nws);
iswin=false(N,nws);
for i=1:nws
    win(:,i)=circshift(win1,startwindow(i)-1);
    iswin(startwindow(i)+(0:nwindow-1),i)=true;
end
c=w'*win./sum(win);
a=(repmat(w,[1 9])-repmat(c,[N 1])).*win; % windowed space/time series, has zero mean

variance=sum(a(:).^2)/N; % windowed variance of time/space series, 88% of variance of k-5/3 spectrum preserved
% Variance is not preserved by the spectral estimation methods, because of
% truncation of the wavenumber space by the Nyquist wavenumber, and
% especially rolloff of the low-wavenumber energy containing range.
% Nevertheless, energy spectral density (per wavenumber) is conserved
% in the middle of the wavenumber range, e.g. knyq/20 - knyq/3.
