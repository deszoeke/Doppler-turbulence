function [sthr,meannoise,kthr]=HildebrandSekhon2(S,R_2cond)
% [sthr,meannoise,kthr]=HildebrandSehkon(F,S)
% White noise threshold and mean finder modified from (Hildebrand and Sekhon 1974) 
% for finding noise beside inertial spectra. Thresholding is along
% frequency dimension, rather than spectral power. R2_cond=2 recommended.
% 2018-05-30 Simon de Szoeke

% Basic idea: exclude large spectral estimates until the remainder
% cannot be distinguished from white noise.
% My method includes lower and lower frequency spectral estimates until one can be
% distinguished from white noise.

ns=length(S)-1; % ignore Nyquist spectral estimate
%[s,ord]=sort(S(1:ns)); % sort s in order of power
ord=(ns:-1:1)';  % order from high frequency to low frequency
s=S(ord);
reord(ord)=1:ns; % indices reverse ordering so s(reord) == S(1:ns)
reord=reord';

if 0<R_2cond<1
% t values for testing whether next point exceeds probability
tt=tinv(R_2cond,(1:length(S))-1); 
% t = (S.-P)./Q % t-like statistic for testing if a single sample belongs
% to the diagnosed distrubtion of white noise.

% Divide each sum by the length of the sum.
% but this results in too conservative a condition, throwing out too many
P=cumsum(s)./(1:ns)'; % <s> for increasingly larger ensembles of s for lower frequencies
% P2=P.*P;
% Q=cumsum(s2)./(1:ns)' - P2;
% R_2=P2./Q;
Q=cumsum((s-P).^2)./(1:ns)'; % <s'^2> for increasingly larger ensembles of larger/low frequency s
R_2=(P.*P)./Q;
% Standard noise threshold is when R_2 goes to 1.
% Collections of spectral estimates S with R_2>1 are white noise;
% collections with R_2<1 are probably signal.
% s(1:ithr) are small noise
ithr=find(R_2<R_2cond,1,'first'); % least nonnoise spectral estimate
if isempty(ithr)
    ithr=ns; % the spectrum is all white noise
end
sthr=s(ithr);

kthr=reord(ithr);    % conservative to exclude kthr from inertial spectrum

% average noise
meannoise=mednmean(s(1:ithr-4),5); % mean of 5 middle points, insensitive to outliers, smoother than median

% Note: HB1974's white noise variance of *frequency* is sigma_N^2 = F^2/12, where
% F=Ndf is the frequency span resolved by the spectrum. This variance
% follows from the spectral density-weighted integral of f'^2 for a flat
% distribution for [0 F]. The span F grows as we include more specral
% estimates. Here we use a nondimensional set of frequencies f.

% commented out to save time, this calculation helps justify the R_1, R_2
% threshold
%{
f=(1:ns)'; % in order
sigN2=f.*f/12; % array increasing in frequency span
css=cumsum(s);
sig2=cumsum(s.*f(ord).*f(ord))./css - (cumsum(s.*f(ord))./css).^2; % HS1974 eqn 4
R_1=sigN2./sig2; % eqn 8

% s2=s.*s;
%}

% R_2cond=1 is too conservative for discerning white
% noise from k^-5/3 turbulence. It cuts off too much turbulence, and sets
% the mean noise too high.

%{
% like HS1974 Fig 2b, but with x-axis reversed
clf
semilogx(s/sthr,R_1,'.-')
semilogx(s/sthr,R_2,'.-')
hold on
plot(1,R_2(ithr),'x')
plot(s([1 end])/sthr,[1 1])
text(0.5,1,'noise')
text(2,1,'signal -->')

clf
loglog(F,S,'.-')
hold on
plot(F(16),S(16),'x')
plot(F([1 end]),meannoise*[1 1])
plot(F([1 end]),sthr*[1 1])

factr=1.3889; F53=F.^(5/3);
vls=factr*(2*pi/Utop)^(2/3)*(F53.*(S-meannoise)); % dissipation ^2/3
epsilon=vls.^1.5; % dissipation
plot(F,epsilon)
plot(F([4 16]),epsilon([4 16]),'o')
plot(F(5:kthr),epsilon(5:kthr),'linewidth',2)
%}

% It looks like the spectrum rolls off at f>F(16)=0.4375 Hz
