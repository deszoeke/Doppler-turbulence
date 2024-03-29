function res = calc_intEps_L(decTime,ts,intNetSpd,USING_TIMESERIES,FETCH_PER_POINT,FETCH_RANGEGATE_RATIO,NUM_LAGS_FIT,specNoiseTailRng,dF,atmVarResult)
% % intEps : dissipation rate (variance method, O'Connor 2010)
% % sigmaW : noise corrected variance
% % L : length scale over N samples (as used in O'Connor)
% % L1 : length scale for individual sample (as used in O'Connor) ; these two describe the scattering volume dimension for the dwell time of the lidar
% % err : noise variance determined from ACF fitting (Lenschow 2000)
% % L0 : found that somewhere without reference (This gets used for Le)
% % Le : instrumental error for integral scale (Lenschow 2000)
% % epsE: instrumental error for dissipation rate (Lenschow 2000, not sure if this is valid (I guess it should be for the variance), but I'm calculating the eps error in a different way (error propagation)
% % L0v : integral length scale from integrating cov up to first zero crossing (f0x)
% intialize the output structure first
%
% Updates added by BJC, Sept 2023, for mobile lidar vertical velocity 
% statistics. Accounts for evenly-spaced data series in space (platform-relative
% fetch) rather than time. I also added the function's input arguments 4 thru [end].
% When using spatial data series, "decTime" is actually an array of distances in
% meters, not time!
% THIS VERSION IS SLIGTHLY MODIFIED FOR DC_FLUX! The fitNoiseXcovSpec_TB(...); uses the arguments from the velStats code that is being used operationally. Also, the version of that fitNoiseXcovSpec_TB() code is likely a bit different than what BJC worked up for mobile vertical velocity statistics
if ~exist('USING_TIMESERIES','var') USING_TIMESERIES = 1; end % Default to the heritage/standard timeseries approach. Alternative is spatial fetch-averaged approach, which was developed for mobile platforms.
if isempty(USING_TIMESERIES) USING_TIMESERIES=1; end
if ~exist('specNoiseTailRng','var') specNoiseTailRng=0.9; end % In case there's an old code that calls this without specNoiseTailRng arg, use the old default value
if isempty(specNoiseTailRng) specNoiseTailRng=0.9; end
if ~exist('atmVarResult','var') atmVarResult=NaN; end
if isempty(atmVarResult) atmVarResult=NaN; end

res.intEps = nan; res.sigmaW = nan; res.L = nan; res.L1 = nan; res.acfNoiseVar = nan; res.L0 = nan; res.Le = nan; res.epsE = nan;
res.L0v = nan; res.acfNorm = nan; res.f0x = nan;

integralScaleIsTime = 0;   % we don't have valid velocity data - express ISL in time
if isnan(intNetSpd) | intNetSpd==0  
    integralScaleIsTime = 1; 
    intNetSpd=1;
end
res.integralScaleIsTime = integralScaleIsTime;
del_spd = 0.25;  % estimating a statistical windspeed uncertainty of 1/4 m/s.

%interpolate TS onto a constant time grid
gdDex = find(~isnan(ts));
if size(gdDex,1) > 1
     ts = ts(gdDex);
     decTime = decTime(gdDex);
     if USING_TIMESERIES  
         dTime = nanmedian(nanOutliers(diff(decTime),1,5));
         numTSteps = (decTime(end)-decTime(1)) / dTime;
         gridTime = decTime(1):dTime:numTSteps*dTime+decTime(1);
         Fs = 1/(dTime*3600*24);  % sample frequency, in 1/seconds
     elseif ~USING_TIMESERIES
         dTime = FETCH_PER_POINT; % Remember, this is not time; it's distance
         numTSteps = (decTime(end)-decTime(1)) / dTime;
         gridTime = decTime(1):dTime:numTSteps*dTime+decTime(1);
         Fs = 1/dTime;  % sample frequency, in 1/meters
     end  
     gridTS =  interp1(decTime,ts,gridTime);
     [xCov,pwrSpec,f] = xCovSpecTB(ts,1,Fs);  % calculate XCOV and spectrum

     % find the inertial subrange - use that to determine range for fit
     first0Dex = find(xCov<0,1,'first')-1;
     if isempty(first0Dex) || first0Dex >= length(f) || first0Dex < 1; return; end
     xCovInertTime = 2.68/(f(end-first0Dex)); % from Lothon 2006, Eq 10. But it's not obviously the same...We don't use L0 so I'll leave this as-is and ignore it! -BJC
     res.L0 = xCovInertTime*intNetSpd/.7468343;  % I think this is some form of Lothon 2006, Eq. 10. L0 is integral scale estimate in the alongwind direction (I guess this refers to single z bin, atmo passing thru horizontally via horz wind. Is the 0.7468 some sort of correction factor? -BJC
     
     startLag = 1;  % Sunil catches this   (was first0Dex);
     stopLag = ceil(Fs*xCovInertTime);
     if startLag >= stopLag; startLag = 1; end
     if stopLag > size(xCov,2)
         stopLag = min(startLag+5,size(xCov,2));
     end
     % startLag is fixed for spatial series, and I've decided to go with a short fixed lag range for fit, as in the run_calcMobileVertVelStats code
     if ~USING_TIMESERIES  
         startLag = FETCH_RANGEGATE_RATIO;
         stopLag = startLag + NUM_LAGS_FIT - 1; %(-1 for indexing to get proper number of lags fit)
     end
     nLagsToFit = stopLag-startLag;
     [varResult] = ...
          fitNoiseXcovSpec_TB(xCov,gridTime,gridTS,Fs,pwrSpec,nanmax([floor(5/(dTime*24*3600)) 3]),5,'',0,0.80,[],dF);

     if varResult.acfNoiseVar < 0 || isnan(varResult.acfNoiseVar)
         varResult.acfNoiseVar = 0;
     end
     
     % Length scale observed per point and per data series (both are rough estimates regarding wind speed)
     varWin = length(gridTime);  % varWin seems to be an artifact from a time when the code might've been calculating var for multiple data windows. As it is now, all the data passed in gets used so varWin is the total data length.
     if USING_TIMESERIES
         T1 = 1/Fs;
         res.L1 = intNetSpd*T1;
         res.L  = intNetSpd*T1*varWin;
     elseif ~USING_TIMESERIES
         res.L1 = FETCH_PER_POINT;
         res.L = FETCH_PER_POINT*varWin;
     end

     % If a variance result was passed in as arg (atmVarResult), use that.
     if ~isnan(atmVarResult)
        res.sigmaW = atmVarResult;
     else
         % Use simple Matlab var() approach to get atmospheric variance estimate
         numBlocks = max(floor(length(gridTime) / varWin),1);
         tVar = []; tVar(1:numBlocks) = NaN;
         for lm = 1:numBlocks
             tts = ts(round(1+(lm-1)*varWin):min(round(lm*varWin),length(ts)));
             tVar(lm) = nanvar(tts);
         end
%          res.sigmaW = tVar-varResult.acfNoiseVar;  % Atmospheric variance
         res.sigmaW = tVar-varResult.lag1Lag0Var;  % Atmospheric variance. Using Lag0-Lag1 because with our modern super-low-instr-noise, extrapolating the ACF does not always work well (e.g. fit projection > lag0). Lag0-Lag1 is guaranteed to behave and give a positive number. This is what Chris tends to use nowadays. -BJC Feb 2024
         if res.sigmaW<0; res.sigmaW=NaN; end
     end
     
     % Calculate integral length scale
     if first0Dex > 1
         res.f0x = first0Dex;
         if USING_TIMESERIES
             x = (0:1/Fs:(first0Dex-1)*1/Fs)'.*intNetSpd; % Array of ACF sample length scale x-axis (in meters), from 0 to ACF zero crossing, in steps of sample period (1/Fs) times the atmosphere-platform speed
         else
             x = (0:1/Fs:(first0Dex-1)*1/Fs)'; % Array of ACF sample length scale x-axis (in meters), from 0 to ACF zero crossing, in steps of sample period (1/Fs)
         end
         % xi, xCovi, ff, fun, y, and zx are all about a 100-point interpolation/fit to the ACF. The values do not change at the original ACF x-points, and I'm not sure the point of this. Maybe it's to interpolate to the zero crossing more precisely? And/or maybe the Matlab "integral" fn does better with more points?
         xi = linspace(x(1),x(end),100)'; % Like x, but 100 evenly-spaced steps, to create a smooth(er) function output later
         xCovi = interp1(x,xCov(1:first0Dex)',xi); % Interpolate xCov to 100 evenly spaced points (xi) between lag 0 and first 0-crossing
         % Calculate the normalized, noise-corrected ACF. BJC is not sure why we do this with the [xCovi(1)-acfNoise] determination of wVar instead of the res.sigmaW we've already calculated, but I'll leave it as-is since it does have a nice symmetry with the xCov noise removal in the numerator, and with sigmaW is may be calc'd from PSD rather than xCov.
         ff = fit(xi,((xCovi-varResult.acfNoiseVar)./(xCovi(1)-varResult.acfNoiseVar)),'spline'); % Fit a spline to the normalized, noise-corrected ACF sampled at xi. THIS NOISE SUBTRACTION MOVES THE ZERO-CROSSING UP!
         fun = @(x) ff(x); % 'fun' is the spline function fit to ACF sampled at the original ACF points
         y = fun(xi); % This is identical to the argument of ff, as far as I can tell. Plotting (xi,y) and (xi,((xCovi-varR.../...NoiseVar))) shows the same thing.
         zx = find(y < 0,1,'first')-1; % zx is index just before the zero crossing. It has changed since first0Dex because of the noise subtraction during the normalization step (ff)
         if isempty(zx) 
             zx = length(x); 
             res.L0v = integral(@(x)fun(x),x(1),x(zx),'ArrayValued',true); % IF there is no zero-crossing in the normalized noise-corrected ACF, take integral of the whole series, which geos to the zero-crossing of the non-normalized, non-noise-corrected ACF.
         elseif zx == 1
             res.L0v = x(2)/2; % IF the normalized noise-corrected ACF is all negative, take the length scale to be half of the lag 1 x-scale (distance period). -BJC: This seems super short! I wonder why they chose this???
         elseif zx > 1
            res.L0v = integral(@(x)fun(x),xi(1),xi(zx),'ArrayValued',true); % Integrate the normalized, noise-corrected ACF interp'd to 100 x-axis points, from lag 0 to the first zero-crossing
         end
     else
         if USING_TIMESERIES
             res.L0v = 1/Fs/2*intNetSpd;
         else
             res.L0v = 1/Fs/2;
         end
     end

     res.intEps = 2*pi*(2/3/0.52)^(3/2)*res.sigmaW^(3/2)*(res.L^(2/3)-res.L1^(2/3))^(-3/2); % O'Connor et al 2010, Eq. 6
     % Regarding above eq': 0.52 is the Kolmogorov constant, which from a quick search of a few references may range 0.50 - 0.55
     
     res.acfNorm = ((xCov-varResult.acfNoiseVar)./(xCov(1)-varResult.acfNoiseVar))'; % was formerly called res.cov. Changed by BJC to be more descriptive; but I don't see this variable getting used anywhere in velStats anyways.
     res.acfNoiseVar = varResult.acfNoiseVar;  % Was formerly res.err. Changed by BJC to be more descriptive; but I don't see this variable getting used anywhere in velStats anyways.

     % I'm not sure how much I trust most of these error fields, since they're downstream of so many of the derived values. And some estimates in there may make a big deal, like del_spd, though I've not sone sensitivity tests.
     res.Le = res.L0v*sqrt(4/size(ts,1)*(res.acfNoiseVar/res.sigmaW));  % lenschow et al 2000, Eq B5. THIS USED TO USE L0 FOR INTEGRAL SCALE, NOT L0v! But I don't understand where L0 is coming from, so I made it L0v. -BJC
     %      res.epsE = res.intEps*sqrt(4/size(ts,1)*(res.acfNoiseVar/res.sigmaW));
     res.sigmaE = res.sigmaW*sqrt(4/varWin*(res.acfNoiseVar/res.sigmaW)); % error for the w variance
     res.epsE = res.intEps*((3*sqrt(res.sigmaE)/sqrt(res.sigmaW))+del_spd/intNetSpd); % propagation of errors for dissipation rate, O'Connor 2010 Eq 11. The del_spd and intNetSpd here are delL and L with the Nt term cancelled out.
end