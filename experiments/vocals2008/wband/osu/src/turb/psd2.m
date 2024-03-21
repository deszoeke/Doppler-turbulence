function [Pxx, Pxxc, f] = psd2(varargin)
%PSD Power Spectral Density estimate.
%   Pxx = PSD(X,NFFT,Fs,WINDOW) estimates the Power Spectral Density of 
%   signal vector X using Welch's averaged periodogram method.  X is 
%   divided into overlapping sections, each of which is detrended, then
%   windowed by the WINDOW parameter, then zero-padded to length NFFT.  
%   The magnitude squared of the length NFFT DFTs of the sections are 
%   averaged to form Pxx.  Pxx is length NFFT/2+1 for NFFT even, (NFFT+1)/2
%   for NFFT odd, or NFFT if the signal X is complex.  If you specify a 
%   scalar for WINDOW, a Hanning window of that length is used.  Fs is the 
%   sampling frequency which doesn't effect the spectrum estimate but is 
%   used for scaling of plots.
%
%   [Pxx,F] = PSD(X,NFFT,Fs,WINDOW,NOVERLAP) returns a vector of frequen-
%   cies the same size as Pxx at which the PSD is estimated, and overlaps
%   the sections of X by NOVERLAP samples.
%
%   [Pxx, Pxxc, F] = PSD(X,NFFT,Fs,WINDOW,NOVERLAP,P) where P is a scalar
%   between 0 and 1, returns the P*100% confidence interval for Pxx.
%
%   PSD(X,...,DFLAG), where DFLAG can be 'linear', 'mean' or 'none', 
%   specifies a detrending mode for the prewindowed sections of X.
%   DFLAG can take the place of any parameter in the parameter list
%   (besides X) as long as it is last, e.g. PSD(X,'mean');
%   
%   PSD with no output arguments plots the PSD in the current figure window,
%   with confidence intervals if you provide the P parameter.
%
%   The default values for the parameters are NFFT = 256 (or LENGTH(X),
%   whichever is smaller), NOVERLAP = 0, WINDOW = HANNING(NFFT), Fs = 2, 
%   P = .95, and DFLAG = 'none'.  You can obtain a default parameter by 
%   leaving it off or inserting an empty matrix [], e.g. PSD(X,[],10000).
%
%   See also CSD, COHERE, TFE, PMEM, PMTM, PMUSIC.
%   ETFE, SPA, and ARX in the Identification Toolbox.

%   Author(s): T. Krauss, 3-26-93
%   Copyright (c) 1988-97 by The MathWorks, Inc.
%   $Revision: 1.24 $  $Date: 1997/02/06 21:56:19 $

%   The units on the power spectra Pxx and Pyy are such that, using
%   Parseval's theorem: 
%
%        SUM(Pxx)/LENGTH(Pxx) = SUM(X.^2)/LENGTH(X) = COV(X)
%
%   The RMS value of the signal is the square root of this.
%   If the input signal is in Volts as a function of time, then
%   the units on Pxx are Volts^2*seconds = Volt^2/Hz.
%
%   Here are the covariance, RMS, and spectral amplitude values of
%   some common functions:
%         Function   Cov=SUM(Pxx)/LENGTH(Pxx)   RMS        Pxx
%         a*sin(w*t)        a^2/2            a/sqrt(2)   a^2*LENGTH(Pxx)/4
%Normal:  a*rand(t)         a^2              a           a^2
%Uniform: a*rand(t)         a^2/12           a/sqrt(12)  a^2/12
%   
%   For example, a pure sine wave with amplitude A has an RMS value
%   of A/sqrt(2), so A = SQRT(2*SUM(Pxx)/LENGTH(Pxx)).
%
%   See Page 556, A.V. Oppenheim and R.W. Schafer, Digital Signal
%   Processing, Prentice-Hall, 1975.

error(nargchk(1,7,nargin))
x = varargin{1};
[msg,nfft,Fs,window,noverlap,p,dflag]=psdchk(varargin(2:end),x);
error(msg)

% compute PSD

window = window(:);
n = length(x);		         % Number of data points
nwind = length(window);    % length of window
if n < nwind               % zero-pad x if it has length less than the window length
   x(nwind)=0;  
   n=nwind;
   disp('The number of points in the time series is less than the window size!');
   disp('Are you sure you want to do this?');
end
x = x(:);	         	      % Make sure x is a column vector
                              % this line must be AFTER the zero-padding in
                              % case x is a scalar
k = fix((n-noverlap)/(nwind-noverlap));     % Number of windows
					               % (k = fix(n/nwind) for noverlap=0)
index = 1:nwind;
KMU = k*norm(window)^2;	% Normalizing scale factor ==> asymptotically unbiased

Spec = zeros(nwind,1); 
Spec2 = zeros(nwind,1);
for i=1:k
   if strcmp(dflag,'none')
        xw = window.*(x(index));
    elseif strcmp(dflag,'linear')
        xw = window.*detrend(x(index));
    else
        xw = window.*detrend(x(index),0);
    end
   index = index + (nwind - noverlap);
    
%*******************************************
%  Changed by JBE  4/3/97
%*******************************************
    xxc=fft(xw,nfft);
    Xx = xxc.*conj(xxc)/nfft;
%    Xx = abs(fft(xw,nfft)).^2;
    Spec = Spec + Xx;
    Spec2 = Spec2 + abs(Xx).^2;
end
  
%*******************************************
%  Added by JBE  4/3/97
%*******************************************

df = Fs/(nfft-1);

% Select first half
if ~any(any(imag(x)~=0)),   % if x is not complex
    if rem(nfft,2),    % nfft odd
        select = (1:(nfft+1)/2)';
    else
        select = (1:nfft/2+1)';
    end
%*******************************************
%  Changed by JBE  4/3/97
%*******************************************
    Spec = 2*Spec(select)/df;
    Spec2 = 4*Spec2(select)/df;
else
    select = (1:nfft)';
%*******************************************
%  Added by JBE  4/3/97
%*******************************************
    Spec = Spec(select)/df;
    Spec2 = Spec2(select)/df;
end
freq_vector = (select - 1)*Fs/nfft;

% find confidence interval if needed
if (nargout == 3)|((nargout == 0)&~isempty(p)),
    if isempty(p),
        p = .95;    % default
    end
    % Confidence interval from Kay, p. 76, eqn 4.16:
    % (first column is lower edge of conf int., 2nd col is upper edge)
    confid = Spec*chi2conf(p,k)/KMU;

    if noverlap > 0
        disp('Warning: confidence intervals inaccurate for NOVERLAP > 0.')
    end
end

Spec = Spec*(1/KMU);   % normalize

% set up output parameters
if (nargout == 3),
   Pxx = Spec;
   Pxxc = confid;
   f = freq_vector;
elseif (nargout == 2),
   Pxx = Spec;
   Pxxc = freq_vector;
elseif (nargout == 1),
   Pxx = Spec;
elseif (nargout == 0),
   if ~isempty(p),
       P = [Spec confid];
   else
       P = Spec;
   end
   newplot;
   plot(freq_vector,10*log10(abs(P))), grid on
   xlabel('Frequency'), ylabel('Power Spectrum Magnitude (dB)');
end
