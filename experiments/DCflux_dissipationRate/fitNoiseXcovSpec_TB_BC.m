function [varResult] = fitNoiseXcovSpec_TB_BC(xcArr,tGrid,tsArr,Fs,psArr,startLag,stopLag,fitOrder,titleStr,plotIt,psNoiseFloorFqCutoff,oneSideFq,dF)
% Brian Carroll updated version of Tim Bonin's version of this code. The code 
% uses the input data timeseries and power spectrum to report various methods
% of variance and noise
% BC update adds TKE dissipation rate (epsilon) to the retrieved values. This is 
% easy to do because the fit to ACF is already being done. Wulfmeyer et al 2024
% has a good discussion of this approach.
% This code is mostly used in calc velStats codes, and in BC's 
% sonicDissipationRate code.
%
% Inputs:  THESE DIFFER FROM TB VERSION b/c I added startLag and stopLag, 
%           instead of just nLagsToFit
% xcArr - autocovariance of the timeseries
% tGrid - timestamps of data series in datenum. Only used in plotting.
% tsArr - data timeseries. Used to plot and to check if multiple, then avg xCov
% Fs - Single value, frequency of data sampling (i.e., inverse of tGrid spacing)
% psArr - Power spectrum of the data series (tsArr)
% startLag - First lag to consider in ACF fit. This is lag number not Matlab
%            array index; i.e., lag 0 is 0 (not 1, tho Matlab starts w/ 1)
% stopLag - Last lag to consider in ACF fit. This is lag number not Matlab ind.
% fitOrder - Nowadays this is just choosing ACFnoise to be fit projection (5) or
%            lag0 - lag1 (0). Used to indicate polynomial fit order. Current 
%            preference is lag0-lag1 b/c fit projection can return negative, but
%            end up reporting 0-1 separately so might as well stick with 5 here.
% titleStr - Can pass in title string for plots, or just be empty as ''
% plotIt - Do plots (1) or no (0)
% psNoiseFloorFqCutoff - Cutoff determines the last percentage of spectrum to 
%                        use for noise floor calculation. The section used is 
%                        1-this, so 0.8 is last 20% for noise floor. Used to be 
%                        named "fc".
% oneSideFq - freq 'x-axis' for power spectrum
% dF - Frequency spacing of the power spectrum (one value, not the freq array)
%
% Edited by BC Feb 2024.


% psNoiseFloorFqCutoff = 0.8;			% cutoff determines last percentage of spectrum to use
% numModes=10;
lngth = length(xcArr(1,:))-1;
if isempty(tsArr) || lngth < stopLag 
    varResult.totSpecVar = nan;
    varResult.totTsVar = nan;
    varResult.totXcVar = nan; 
    varResult.specNoiseVar = nan;
    varResult.acfNoiseVar = nan;
    varResult.lag1Lag0Var = nan;
    varResult.numTS = nan;
    varResult.numPts = 0;
    varResult.ACF_fitCoef = nan;
    varResult.ACF_fitCoef_confBounds = nan;
%    varResult.GWvar=nan;
    
%     GW.wvar(1:numModes,1)=NaN;
%     GW.avg_amp(1:numModes,1)=NaN;
%     GW.avg_period(1:numModes,1)=NaN;
  return
end
 
% for plotting
f = Fs/2*linspace(0,1,length(psArr(1,:))); % x axis

%% If we have more than one time series in the arrays, average them together.
multiArrays = min(size(tsArr)) > 1;  
if multiArrays
    mnXcArr = mean(xcArr);
    meanPwrSpec = mean(psArr);
else
    mnXcArr = xcArr;
    meanPwrSpec = psArr;
end

%% Old GW detection stuff. Was already commented out before BC interacted with the code.
% First need to remove gravity waves and quantify their impact
% xCov=mnXcArr;
% avgDifft=1./Fs;
% if avgDifft>10
%     avgDifft=10;
% end
% minPtsForGWDetection=round(60./avgDifft); %need at least 1 minutes of data
%numModes=5;
% gofThresh=0.95; % Needs to be higher than this to be a gravity wave
% GW.wvar(1:numModes,1)=NaN;
% GW.avg_amp(1:numModes,1)=NaN;
% GW.avg_period(1:numModes,1)=NaN;
% isGWPresent=1;
% 
% x1Temp=find(oneSideFq>0.1);
% HFvar=nansum(meanPwrSpec(x1Temp)).*dF;
% x1Temp=find(oneSideFq<=0.1 & oneSideFq>0.01667);
% MFvar=nansum(meanPwrSpec(x1Temp)).*dF;
% x1Temp=find(oneSideFq<=0.01667);
% LFvar=nansum(meanPwrSpec(x1Temp)).*dF;
% [PK2,~]=findpeaks(smooth(xCov));
% [LP1,~]=findpeaks(-smooth(xCov));
% ALLFVar=LFvar+MFvar+HFvar;
% if~isempty(LP1)
%     LP1=LP1(1);
% else
%     LP1=0;
% end
% if ~isempty(PK2)
%     PK2=PK2(1);
% else
%     PK2=0;
% end
% 
% if HFvar+MFvar<0.02 & LFvar>0.01
%     isGWPresent=1;
% %elseif HFvar+MFvar<0.2*ALLFVar & MFvar+HFvar<0.1 %hard limit
% %    isGWPresent=1;
% else
%     isGWPresent=0;
% end
% 
% if length(xCov)>minPtsForGWDetection
% for i=1:numModes % first 5 modes.... should remove nearly all variance from GWs
%     if   xCov(1)/xCov(2)>2 || xCov(1)<0.005 % || diff(xCov(2:3))<diff(xCov(firstZeroCrossing-1:firstZeroCrossing))%|| firstZeroCrossing/2<halfCrossing %nanmean((diff(xCov(round(firstZeroCrossing/2):firstZeroCrossing))))>nanmean((diff(xCov(1:round(firstZeroCrossing/2))))) 
%        isGWPresent=0; 
%     end
% 
%     
%     if isGWPresent
%     firstZC=find(xCov<0,1);
%     if isempty(firstZC)
%         [~,firstZC]=nanmin(xCov);
%     end
%     if firstZC>4 %skip if it is small... not a GW
%     %initialIndxToUse=[ceil(minPtsForGWDetection/20):floor(minPtsForGWDetection/4)];
%     initialIndxToUse=[1:firstZC];
% if nansum(GW.wvar)==0
% s = fitoptions('Method','NonlinearLeastSquares',...
%     'Lower',[0 0 0],...
%     'Upper',[10 5 0.0001],...
%     'Startpoint',[xCov(2) pi/2/firstZC 0]);    
% else
% s = fitoptions('Method','NonlinearLeastSquares',...
%     'Lower',[0 0 -1],...
%     'Upper',[10 5 1],...
%     'Startpoint',[xCov(2) pi/2/firstZC -nansum(GW.wvar)]);    
% end
% ft=fittype('a*cos(x*b)+c','options',s);
% [ACFfit,gof2]=fit((initialIndxToUse-1)',xCov(initialIndxToUse)',ft);
% resid=abs(xCov-(ACFfit.a*cos(ACFfit.b*[0:length(xCov)-1])+ACFfit.c));
% %refinedIndx=find(resid(1:minPtsForGWDetection)<nanmedian(resid(1:minPtsForGWDetection)));
% refinedIndx=find(resid(1:minPtsForGWDetection)<nanmedian(resid(1:minPtsForGWDetection)));
% refinedIndx=refinedIndx(1:nanmin([20 length(refinedIndx)])); %only use 20 points for fit
% [ACFfit,gof2]=fit((refinedIndx-1)',xCov(refinedIndx)',ft);
% % s = fitoptions('Method','NonlinearLeastSquares',...
% %                'Lower',[],...
% %                'Upper',[],...
% %                               'Startpoint',[1 0]);
% % %                'Lower',[0,0],...
% % %                'Upper',[Inf,max(stopLag)],...
% % 
% % % % % ft = fittype('a*(x-b)^n','problem','n','options',s);
% %  ft = fittype('a*(x)^n+b','problem','n','options',s);
% %  [~,gof3] = fit((refinedIndx-1)',xCov(refinedIndx)',ft,'problem',2/3);
% 
% if gof2.rsquare>gofThresh
% GW.wvar(i,1)=ACFfit.a+ACFfit.c;
% GW.avg_amp(i,1)=2*ACFfit.a;
% GW.avg_period(i,1)=2*pi./ACFfit.b./Fs;
% xCov=xCov-(ACFfit.a*cos(ACFfit.b*[0:length(xCov)-1])+ACFfit.c);
% else
%     isGWPresent=0;
% end
%     end
%     end
% end
% end
% %mnXcArr=xCov;
%  GWvar=nansum(GW.wvar);

%% Fit the ACF with -2/3 power tau dependency, e.g. Wulfmeyer et al. 2024
%stopLag = nanmax([stopLag 3]);
% startLag = 1;

y = mnXcArr(startLag+1:stopLag+1);
x = [startLag:stopLag];
X = [0:stopLag];  %Big X and Y are just for plotting
Y = mnXcArr(1:stopLag+1);

% Fit ACF to -2/3 power law
s = fitoptions('Method','NonlinearLeastSquares',...
               'Lower',[],...
               'Upper',[],...
                              'Startpoint',[1 0]);
%                'Lower',[0,0],...
%                'Upper',[Inf,max(stopLag)],...

 ft = fittype('a*(x)^n+b','problem','n','options',s);
 [ACFfit,gof2] = fit(x',y',ft,'problem',2/3);

 % Y-intercept of the fit is lag0, the atmospheric variance
 A(1)=ACFfit.b;
 if A(1)<0
    A(1)=0.000000001; % Making this very small instead of 0, since 0 can throw some errors in codes downstream
 end

 % ACF fit coefficient is related to the dissipation rate, so we output it.
 ACF_fitCoef = -ACFfit.a;
 temp = -confint(ACFfit);
 ACF_fitCoef_confBounds = temp(:,1);
 % The confidence bounds should be symmetric, I think. So just report the difference as a fit confidence error.
 ACF_fitCoef_cbError = abs(abs(ACF_fitCoef_confBounds(1))-ACF_fitCoef);
 
 % Simple diagnostic plot of ACF, lags selected for fit, and the fit.
% % figure(22);clf;
% % plot(X,Y,'o',x,y,'.');
% % hold on
% % plot(ACFfit,'m')

% Old code - polynomial fit instead of -2/3 power structure function:
% A = polyfit(x,y,fitOrder);

lag1Lag0Var = mnXcArr(1)-mnXcArr(2);    % form difference of 1st 2 pts

%% Get ACF noise variance and X, Y for plotting
switch fitOrder
    case 0  % don't do a fit
        Y = X;  % to get the size right
        Y(:) = mnXcArr(2);
        acfNoiseVar = lag1Lag0Var;    % form difference of 1st 2 pts
    case 5 %actually structure function fit
        Y= X; %just a filler
        acfNoiseVar=mnXcArr(1)-A(1);
        % Old code - polynomial fit stuff
%     case 1
%         Y = A(1)*X+A(2);
%         acfNoiseVar = mnXcArr(1)-A(2);
%     case 2
%         Y = A(1) * X.^2+ A(2) * X.^1  + A(3);
%         acfNoiseVar = mnXcArr(1)-A(3);
%     case 3
%         Y = A(1) * X.^3+ A(2) * X.^2 + A(3)*X + A(4);
%         acfNoiseVar = mnXcArr(1)-A(4);
%     case 4
%         Y = A(1) * X.^4+ A(2) * X.^3 + A(3)*X.^2+A(4)*X + A(5);
%         acfNoiseVar = mnXcArr(1)-A(5);
end


% Spectral noise
nsStopDex = length(psArr(1,:));
nsStartDex = nsStopDex-floor((1-psNoiseFloorFqCutoff)*nsStopDex);
specNoiseAvg = mean(meanPwrSpec(nsStartDex:nsStopDex));
deltaF = 1/(range(tGrid)*24*3600);
specNoiseVar = specNoiseAvg * nsStopDex * deltaF;    

% Build the output structure
varResult.totSpecVar = sum(meanPwrSpec)*dF;%deltaF;       % the integral is dF*sum(Pi)
varResult.totTsVar = mean(nanvar(tsArr'));
varResult.totXcVar = mnXcArr(1);
varResult.specNoiseVar = specNoiseVar;
varResult.acfNoiseVar = acfNoiseVar;
varResult.lag1Lag0Var = lag1Lag0Var;
varResult.numTS = length(tsArr(:,1));
varResult.numPts = length(tsArr(1,:));
varResult.ACF_fitCoef = ACF_fitCoef;
varResult.ACF_fitCoef_confBounds = ACF_fitCoef_cbError;
%varResult.GW=GW;
%varResult.GWvar=GWvar;
 
%% Plot results. Multiple subplots showing data series, ACF, power spectrum, etc
meanPwrSpecSmoothed = smooth(meanPwrSpec,3);
if plotIt
    fig = figure(33);
    fig.Position = [16 60 1210 565];

    % Data timesries
    subplot(3,2,[1,2])
    plot((tGrid-tGrid(1))*24*3600,tsArr,'.-')
    title(titleStr);
    % datetickzoom;
    xlabel('Relative Time (sec)')
    ylabel('Velocity (m/s)');
    grid on

    % Zoomed-out view of ACF and fit
    subplot(3,2,3)
    plot([0:lngth],xcArr)
    if multiArrays
        hold on
        plot([0:lngth],mnXcArr,'ko' ,'MarkerFaceColor','k');
        hold off
    end
    hold on
    plot(ACFfit,'r')
    grid on
    xlim([0 100])% 300   %length(xcArr)])   % stopLag*10
    ylim([-0.01 inf])
    xlabel('Lag #')

    ylabel('ACF ');
    title(['ACF fit to lags ' num2str(startLag) ' to ' num2str(stopLag) '     Measurement freq: ' num2str(Fs,3) ' Hz'])
    hold off

    % Zoomed-in view of the ACF and fit from lag0 to 2*stopLag
    subplot(3,2,5)
    hold off
    plot([0:2*stopLag],mnXcArr(1:2*stopLag+1),'kd');
    hold on
    plot([startLag stopLag],[mnXcArr(startLag+1) mnXcArr(stopLag+1)],'*','Color','blue');
    xlabel('Lag #')
    ylabel('ACF ');
%     plot(X,Y)
    plot(ACFfit,'r')
    xlim([0 2*stopLag])
    ylim([nanmin(mnXcArr(1:2*stopLag+1)) mnXcArr(1)+0.01])
%     titleStr = sprintf('ACF Noise stdev : %1.2f    Slope : %1.2f/3',sqrt(acfNoiseVar),A(1)*3);
    titleStr = sprintf('ACF Noise Var: %1.3f      lag0-lag1 Var: %1.3f',acfNoiseVar,lag1Lag0Var);
    title(titleStr);
    grid on
    hold off

    % Power spectrum and such
    subplot(3,2,[4,6])
    hold off
    loglog(f,meanPwrSpecSmoothed,'.');
%         loglog(f,meanPwrSpec,'.');
    hold on
    loglog([f(2),f(end)],[specNoiseAvg,specNoiseAvg],'--g');
    loglog([f(nsStartDex),f(nsStartDex)],[min(meanPwrSpec(2:end)),max(meanPwrSpec(2:end))],'--');
    loglog([f(nsStopDex),f(nsStopDex)],[min(meanPwrSpec(2:end)),max(meanPwrSpec(2:end))],'--');
    % loglog([fFit1,fFit1],[min(exp(y)),max(exp(y))],'-.k');
    % loglog([fFit2,fFit2],[min(exp(y)),max(exp(y))],'-.k');
    % loglog(exp(x),exp(y),'r-.');
%     titleStr = sprintf('Spectral Noise stdev : %1.3f    Slope : %1.2f/3',sqrt(specNoiseVar),numerator);
    titleStr = sprintf('Spectral Noise Var : %1.3f ',(specNoiseVar));
    title(titleStr);
    xlabel('Frequency (Hz)')
    ylabel('Power (arb units)');
    xlim([nanmin(f(2:end)) nanmax(f(2:end))])
    ylim([nanmin(meanPwrSpecSmoothed(2:end)) nanmax(meanPwrSpecSmoothed(2:end))])

    % -5/3 power reference line
    referenceScaling = specNoiseAvg * f(nsStartDex)^(5/3); % plotting the -5/3 line uses y=a*x^-5/3, where a scales to be near data magnitude. Find a by using a known point of interest, which in this case is where the noise floor begins.
    referenceX = f; %f(floor(length(f)/20):end);  % latter 95% of the x-points (frequencies)
    referenceY = referenceScaling .* referenceX.^(-5/3);
    loglog(referenceX,referenceY,'r--')

    subplot(3,2,[1,2]) % Make this one the active axes so we can set title to it in run_calcMobileVertVelStats
%     title(['Hoover building sonic anemometer   ' datestr(tGrid(1),'yyyy-mm-dd HH:MM:SS') '  -  ' datestr(tGrid(end),'HH:MM:SS')])
%     print(['C:\Users\bcarroll\Documents\2022_DCFLUX\UrbanNet_HooverTowers\images\ACFplots\HooverBldgSonicACFplot_' datestr(tGrid(1),'yyyymmdd_HHMMSS')],'-dpng')
    % Need to print lidar plots outside of the loop so we can have height info.
%     title(['Hoover building lidar   ' datestr(tGrid(1),'yyyy-mm-dd HH:MM:SS') '  -  ' datestr(tGrid(end),'HH:MM:SS')])
%     print(['C:\Users\bcarroll\Documents\2022_DCFLUX\halo49\images\ACFplots\HooverBldgLidarACFplot_' datestr(tGrid(1),'yyyymmdd_HHMMSS')],'-dpng')
end