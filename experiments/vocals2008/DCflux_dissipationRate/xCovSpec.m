function  [xCov,pwrSpecDensity,oneSideFreq] = xCovSpec(ts,rmvMean,Fs)
% calculates crosscovariance and power spectrum of time series
if ~exist('rmvMean','var') rmvMean = 1; end   % assume we want to unless specifically told no ...
if ~exist('Fs','var') Fs = 1; end   
if rmvMean ts = ts - (nanmean(ts)); end
nPts = length(ts);
% The covariance
tXcov = xcov(ts);%/nPts;

% The spectral calculation
ftTs = fftshift(fft(ts,nPts))/nPts; % calc the ft
pSpec = abs(ftTs.*conj(ftTs)); % two sided powerspectrum
f=[-nPts/2+1:nPts/2]/(nPts/Fs);
% One sided version
cntrPt = floor(nPts/2) + 1;  %position of the DC offset in the shifted spectrum
pwrSpec = pSpec(cntrPt:end);  % transfer the dc value (cntrPt) and the ultimate size of the array
if mod(nPts,2) lastPoint = 1; else lastPoint = 2; end  % if even vs odd
pwrSpec(2:end)= pSpec(cntrPt+1:end) + pSpec(cntrPt-1:-1:lastPoint);  % single sided power density (sum up pos and neg freqs)
oneSideFreq = f(cntrPt:end);    % this is the actual frequency (not a unit frequency)
pwrSpecDensity = pwrSpec /mean(diff(oneSideFreq));  % to make it a density - divide by dF (normalized frequency)

% now the xcov
cntrPt = floor(length(tXcov)/2)+1;
xCov = tXcov(cntrPt:end);
[m1,n1] = size(xCov);
den = (nPts-[0:length(xCov)-1]);

[m2,n2] = size(den);
if (m1==m2 && n1==n2)
    xCov = xCov./den;
elseif (m1==n1 && m2==n2)
    xCov = xCov./den';
end