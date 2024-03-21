function [sthr,meannoise,kthr]=HildebrandSekhon(S)
% [sthr,meannoise,kthr]=HildebrandSehkon(F,S)
% Hildebrand and Sekhon (1974) white noise threshold and mean.
% 2018-05-30 Simon de Szoeke

% Basic idea: exclude large spectral estimates until the remainder
% cannot be distinguished from white noise.
% Here we include larger and larger spectral estimates until one can be
% distinguished from white noise.

ns=length(S)-1; % ignore Nyquist spectral estimate
[s,ord]=sort(S(1:ns));

% HB1974's white noise variance of *frequency* is sigma_N^2 = F^2/12, where
% F=Ndf is the frequency span resolved by the spectrum. This variance
% follows from the spectral density-weighted integral of f'^2 for a flat
% distribution for [0 F]. The span F grows as we include more specral
% estimates. Here we use a nondimensional set of frequencies:
f=(1:ns)'; % in order
sigN2=f.*f/12; % array increasing in frequency span
css=cumsum(s);
sig2=cumsum(s.*f(ord).*f(ord))./css - (cumsum(s.*f(ord))./css).^2; % HS1974 eqn 4
R_1=sigN2./sig2; % eqn 8

s2=s.*s;

% Dividing all sums by ns is wrong, but
% it seemed to give good noise estimates
% P=cumsum(s)/ns; 
% Q=cumsum(s2)/ns - P2;
% Corrected: divide each sum by the length of the sum.
% but this results in too conservative a condition, throwing out too many
P=css./(1:ns)'; % <s> for increasingly larger ensembles of larger s
% P2=P.*P;
% Q=cumsum(s2)./(1:ns)' - P2;
% R_2=P2./Q;
Q=cumsum((s-P).^2)./(1:ns)'; % <s'^2> for increasingly larger ensembles of larger s
R_2=(P.*P)./Q;
% Threshold is when R_2 goes to 1.
% Collections of spectral estimates S with R_2>1 are white noise;
% collections with R_2<1 are probably signal.
% s(1:ithr) are small noise
ithr=find(R_2<1,1,'first'); % least nonnoise spectral estimate
if isempty(ithr)
    ithr=ns; % the spectrum is all white noise
end
sthr=s(ithr);

% exclude first S(k) to cross below threshold, and exclude higher frequencies
kthr=find(S<sthr,1,'first')-1;
meannoise=mean(s(1:ithr));

% This is procedure is right, but too conservative for discerning white
% noise from k^-5/3 turbulence. It cuts off too much turbulence, and sets
% the mean noise too high.

% some other metrics might be helpful
mednoise=median(s(1:ithr));
% try cutting the spectrum by frequency, rather than by threshold

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
