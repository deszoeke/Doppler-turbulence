function funcCalcHaloVertVelStats_epsilon(params,workDate,addOn,totHaloData)
% --------------------------------------------------------------------------
% Calculates vertical velocity stats
% pulls in boundary layer heights and Vert beam info - use this to display
% vert data
%
% 4/28/13 WAB
%
% 3/07/15 TAB - Updated to calculate noise through ACF method, and quantiy
% HF turbulence used for BLH
% 
% March 2024 - Brian Carroll - Updated to use ACF fit coefficient for 
% calculation of dissipation rates (epsilon). Whomever did this before me was
% using O'Connor et al 2010 method, which I also brushed up a bit months ago, 
% but consensus is towards the ACF fit coef method now (Wulfmeyer et al 2024)
%   This change is largely within the subfunction fitNoiseXcovSpec_TB_BC
%--------------------------------------------------------------------------

% workDate= datenum([2013,04,25,00,0,0]);
% params = influxParams;

if nargin<4 %if totHaloData not defined
    totHaloData=[];   % BJC - I have no idea what this is! It's empty in the deployed Halo process routine, so it'll stay empty here. :)
end


warning off
expStr = params.expStr;
expDir = params.expDir;
diskTag = params.diskTag;
scannerHeight = params.scannerHeight;            % how high the scanner is over the water
onlyUseCloudTopped = params.vertVelStatsOnlyUseCloudTopped;  % don't calc for other scenarios
maxTotVarThresh = params.vertVelStatsMaxTotVarThresh; % totvar larger than this - don't bother
minLegitElevation = params.vertVelStatsMinLegitElevation;

outputDataDir = params.velStatsDir;
imageOutputDir = [params.imageDir,'velStats\'];
movieOutputDir = params.movieDir;

HFvarFqThreshold=0.01667; %frequencies higher than this are considered HF

if ~isfield('velStatsTitleStr','params')
    titleStr = 'HALO Vertical Data';
else
    titleStr = params.velStatsTitleStr;
end

version = 120;

dumpGraphics = 1;  % This is inside the IF for displayGraphics, so can just always leave as 1 and use displayGraphics to control plot production and saving.
displayGraphics = 0;
oneMovieFile = 0;
makeMovie = 0;
PLOT_ACF = 0;   ALTERNATE_COUNTER=0; %Leave ALTERNATE at 0; used as an alternating flag to only plot every other altitude (save time and space)
outputData = 1;
removeRain = 0;
rainVelThresh = -2;  % applied to column/row average
CALC_EPSILON_AND_L_OC = 1; % Calculate dissipation rate
loadWindProfiles = 1; % For dissipation rate (epsilon) and  integral length scale
WIND_SPD_UNC = 0.75;  % Estimate of uncertainty in wind speed (m/s). This is probably a bit high for calm conditions, but maybe decent for daytime CBL.
                      % Ideally we'd use unc propoagated from radial vel and/or output from VAD process, but weirdly that doesn't exist for the DCFlux dataset I'm currently working with. Not going down that rabbit hole right now.

% loadSmoothedWindProfs = 0;
workLoader;
snrLims = [-25,5];

mov = [];

yMax = params.vertVelStatPlotMaxY;
numYTicks = 5;
dY = 8; % pixels for small spacing on the bottom of the upperplot in the Y
dX = 70; % pixel spacing on left and right
extraYTopRoom = 15; % extra spacing at top for upper plot
lessYBotRoom = 10; % shave and stretch lower plot by this positive amount
spc    = '                                                       ';
lilSpc = '                          ';
lilSpc = '   -   '; % scaled down for png output
if displayGraphics
    figure(301)
    set(gcf,'units','norm','pos',[0.0805    0.1333    0.8742    0.7750])
    setCarboneColormap;
end

flag = -9999;
[ideal_dCRLB,pCRLB2,crlbSNRs] = calcIdealCRLB(expStr);

outFileStr = 'velStats';

oldMovDay = 0;          % start at zero to trigger a new movie file on first round
oldDecDate = -1;
outRecCount = 0;
clear velStats;

% nanVarResult = fitNoiseXcovSpec_TB([1],[],[],0,[],3,1,'',0,0);  % pull a nan structure out
nanVarResult = fitNoiseXcovSpec_TB_BC([1],[],[],0,[],3,5,1,'',0,0);  % pull a nan structure out
nanEpsilon = calcEpsilon([],[],0,0,0,0);
outFilePending = 0;

qp = params.patternRepeatPeriod;
[haloDataArr,haloData,totHaloData] = selectHaloDataByScanType(workDate,workDate+1,[0],params,totHaloData);
if isempty(haloDataArr) return; end;
haloDataArr = calcHaloFileNums(haloDataArr);

if addOn
   startDate = workDate;
   stopDate = workDate+1;
   loadVelStats = 1;
   workLoader;
   if exist('velStats','var') &  numel(velStats) > 2
       velStats = velStats(1:end-1);  % clip off end and re-calc the last one
       existingVelStatsEndTime = velStats(end).decTime;
   else
       addOn = 0;
   end
else
    clear velStats
end

% Load in wind profiles for later integral scale calculation
if loadWindProfiles
    loadWindProfs = 1;   loadVelStats = 0;    startDate=workDate;  stopDate=workDate+1;
    workLoader;          loadWindProfs = 0;   startDate=[];         stopDate=[];
    % wind profiles get loaded in as windProfs with a struct index for each flt&sp that was loaded. We need to concatenate all of them into one, called windProf, to interface with the rest of the code as-is.
    if ~exist('windProf','var') 
    windProf.decTime=[];   windProf.u=[];   windProf.v=[];   windProf.height=[];
    for i = 1:length(windProfs)
        windProf.decTime = vertcat(windProf.decTime,nanmean(windProfs(i).avgTime));
        windProf.u = horzcat(windProf.u,windProfs(i).u);
        windProf.v = horzcat(windProf.v,windProfs(i).v);
        windProf.height = horzcat(windProf.height,windProfs(i).avgZ);      % Height here is AGL; make sure it agrees with any other heights in this code! (We store MSL too so easy to change if needed)
    end
    windProf.height = median(windProf.height,2,'omitnan');
%     windProf.u=windProf.u'; windProf.v=windProf.v'; windProf.height=windProf.height';%windProf.decTime=windProf.decTime';  %need to transpose to fit with preexisting code
    end
end

% run the times on queue period boundaries
for testTime = floor((workDate)/qp)*qp : qp : ceil((workDate+1)/qp)*qp-qp;
    startTime = testTime;
    stopTime = testTime + qp;
    if addOn & (stopTime < existingVelStatsEndTime)
        fprintf('Adding on to end - skipping %s\n',datestr(stopTime));
        outRecCount = outRecCount + 1;
        continue;
    else
        fprintf('Calculating for %s - %s\n',datestr(startTime,'HH:MM'),datestr(stopTime));
    end
    % HALO's Queue is not synchronized to UTC so find the true imits of the
    % file we want independent of the start/stop time values
    dex = find(haloDataArr.decTime >= startTime & haloDataArr.decTime < stopTime & haloDataArr.elevation>minLegitElevation);
    if isempty(dex) continue; end
    medianFileNum = round(nanmedian(haloDataArr.fileNum(dex)));
    dex = find(haloDataArr.fileNum == medianFileNum);
    if numel(dex) < 5 continue; end
  
    decTime = haloDataArr.decTime(dex);
    ranges = haloDataArr.range;
    vel = haloDataArr.velocity(:,dex);
    snr = haloDataArr.snr(:,dex);
    beta = haloDataArr.beta(:,dex);
    elev = haloDataArr.elevation(dex);
%     startTime = decTime(1);
%     stopTime = decTime(end);
    
    fprintf('Actual :  %s - %s\n',datestr(decTime(1),'HH:MM'),datestr(decTime(end)));
    
    jan0th = calcJan0th(decTime);
    julianDay = floor(startTime-jan0th);
    
    meanLat = params.siteLat;
    meanLon = params.siteLon;
    meanDecTime = nanmean(decTime);

    fprintf('Working Date : %s\n',datestr(meanDecTime));

    % Get a wind profile for integral scale. Taking wind profile interpolated to the mean decTime of the velStats vertical stare data
    if loadWindProfiles
        % Interpolate wind profile timeseries to get a wind profile at the center of this velStats timespan
        clear intWpU intWpV
        for jk = 1:length(windProf.height)
            gdDex = find(~isnan(windProf.decTime') & ~isnan(windProf.u(jk,:)));
            if length(gdDex)>1 % Require at least two wind measurements. Could allow just one if it's close enough in time, but if there's only one than we're likely in a noisy regime so skip it.
                % Make sure one of the wind profiles is within at least 40min of query time, so we're not using data from super far away
                if min(abs(windProf.decTime(gdDex)-meanDecTime)) < datenum(0,0,0,0,40,0)
                    intWpU(jk) = interp1(windProf.decTime(gdDex),windProf.u(jk,gdDex),meanDecTime);
                    intWpV(jk) = interp1(windProf.decTime(gdDex),windProf.v(jk,gdDex),meanDecTime);
                else
                    intWpU(jk) = NaN;
                    intWpV(jk) = NaN;
                end
            else
                intWpU(jk) = NaN;
                intWpV(jk) = NaN;
            end
        end
        % Interpolate the wind profile to vertical stare heights
        intU = nanInterp1(windProf.height,intWpU,ranges,0);
        intV = nanInterp1(windProf.height,intWpV,ranges,0);
        [spdProf,dirProf] = rthetaFromUV(intU,intV);
    end

    badDex = find(elev < minLegitElevation);
    vel(:,badDex) = nan;
    
    % apply SNR Threshold : vertVelStatsMinSnrThresh
    badDex = find(snr < params.vertVelStatsMinSnrThresh);
    vel(badDex) = nan;
    
    
%     % look for rain ...
%     if removeRain
%         % pull out proper column
%         [v,useDex] = min(abs(nanmean(decTime)-meanProf.meanDecTime));
%         % we are using rows and cols - but the Vd array is transposed to
%         % start with ...
%         [r,c] = size(vel);  % these arrays are oversized ...
%         [tR,tC] = size(meanProf.meanVelProf(useDex,:));
%         
%         badRowDex = find(abs(meanProf.meanVelProf(useDex,1:min([c,tC]))) > abs(rainVelThresh) );
%         [tR,tC] = size(meanProf.colMeanVelProf(useDex,:));
%         badColDex = find(abs(meanProf.colMeanVelProf(useDex,1:min([r,tC]))) > abs(rainVelThresh) );
%         % we can also get large positive stretches of data associated
%         % with a hard target in the noise gates - using abs values
%         % should pick that up as well
%         Vd2 = Vd;
%         Vd2(:,badRowDex) = nan;
%         Vd2(badColDex,:) = nan;
%         if sum(size(Vd)) ~= sum(size(Vd2))
%             a=1;
%         end
%         figure(100)
%         plot(decTime,meanProf.colMeanVelProf(useDex,1:min([r,tC])),...
%             decTime(badColDex),meanProf.colMeanVelProf(useDex,badColDex),'rs',...
%             decTime,nanmean(Vd'),'g'); datetick('keeplimits');
%         figure(102)
%         plot(meanProf.meanVelProf(useDex,1:min([c,tC])),ranges,meanProf.meanVelProf(useDex,badRowDex),ranges(badRowDex),'rs',nanmean(Vd),ranges);
%         ylim([0,3000]);
%         Vd = Vd2;
%         
%     end

% transfer the nans ...
    vel(isnan(snr)) = nan;
    snr2=snr; % don't remove nans here... use all values for calculating variance of SNR
    snr(isnan(vel)) = nan;
    beta(isnan(vel)) = nan;
    
    
%     CI = calcCloudInfo(VBI(jk),beamInfo,minLegitElevation);
%     if onlyUseCloudTopped
%         nanDex = ones(nb,1);
%         nanDex(CI.posNegGradPeakDex) = 0;  % going to use this to nan out bad data
%     else
%         nanDex = zeros(nb,1);      % initialize to zero if not using it
%     end
    ppb = numel(ranges);
    
    avgSnr = nan(1,ppb);
    varSnr = nan(1,ppb);
    avgBeta = nan(1,ppb);
    meanResult = nan(1,ppb); % this is the mean with col avg removal (if any)
    noOffsetMeanResult = nan(1,ppb);  %This is to restore the offset from the col avg removal (if any)
    varResult = nan(1,ppb);
    skewResult = nan(1,ppb);
    kurtResult = nan(1,ppb);
    crlbVarResult = nan(1,ppb);
    xCovVarResult = nan(1,ppb);
    pwrSpecVarResult = nan(1,ppb);
    HFvar = nan(1,ppb);
    pctElHFvar = nan(1,ppb);
    nResult = nan(1,ppb);
    epsilon = nan(1,ppb);
    epsilonError = nan(1,ppb);
    epsilon_OC = nan(1,ppb);
    epsilonError_OC = nan(1,ppb);
    integralScaleResult = nan(1,ppb);
    integralScaleError = nan(1,ppb);
    integralScaleIsTime = nan(1,ppb);
    clear gridVarResult

    for kl = 1:ppb
        gdDex = find(~isnan(vel(kl,:))); % calculate across in time
        ts =  vel(kl,gdDex);
        is =  snr(kl,gdDex);
        is2=  snr2(kl,:);
        ib =  beta(kl,gdDex);
        t = decTime(gdDex);
        % remove outliers in velocity ....
        ts = nanOutliers(ts,2,5);  % 3 sigma , 2 passes
        gdDex =  find(~isnan(ts));
        [~,gdDex2]= unique(decTime(gdDex)); %removing duplicates in time (if any)... needed for interpolation later
        gdDex=gdDex(gdDex2);
        ts = ts(gdDex);
        is = is(gdDex);
        ib = ib(gdDex);
        t = t(gdDex);
        avgSnr(kl) = 10*log10(nanmean(10.^(is2/10)));
        varSnr(kl) = nanvar(10.^(is2/10));
        
        if length(ts) > 20 & nanvar(ts) < maxTotVarThresh
   %         avgSnr(kl) = 10*log10(nanmean(10.^(is/10)));
            avgBeta(kl) = 10*log10(nanmean(10.^(ib/10)));
            %interpolate TS onto a constant time grid
            dTime = nanmean(nanOutliers(diff(t),1,3));
            numTSteps = (decTime(end)-decTime(1)) / dTime;
            gridTime = [decTime(1):dTime:numTSteps*dTime+decTime(1)];
            dTimeSec = (gridTime-gridTime(1))*3600*24;
            Fs = 1/(dTime*3600*24);
            gridTS =  interp1(decTime(gdDex),ts,gridTime);
            
            [xCov,pwrSpec,oneSideFq] = xCovSpec(ts,1,Fs);  % calculate XCOV and spectrum
            dF=nanmean(diff(oneSideFq));
            x1Temp=find(oneSideFq>HFvarFqThreshold);
            HFvar(kl)=nansum(pwrSpec(x1Temp)).*dF;
            pctElHFvar(kl) = length(x1Temp)./ length(oneSideFq);
            x1Temp=find(oneSideFq<=HFvarFqThreshold);
       %     LFvar(kl,outRecCount+1)=nansum(pwrSpec(x1Temp)).*dF;
%             gridVarResult(kl) = fitNoiseXcovSpec_TB(xCov,gridTime,gridTS,Fs,pwrSpec,nanmax([floor(5/(dTime*24*3600)) 3]),5,'',0,0.80,oneSideFq,dF);
            % Bonin et al 2016, range gate size divided by horz spd. This lower limit avoids the eddies that are too small to be resolved and would thus be underestimated variance.
            startLag = min([75 round( nanmean(diff(ranges)) ./ spdProf(kl) )]); % min([75 x]) to make sure we don't get some crazy number that breaks the code if winds are near zero. Reevaluate this based on 
            stopLag = startLag + round( nanmax([floor(5/(dTime*24*3600)) 3]) );
            gridVarResult(kl) = fitNoiseXcovSpec_TB_BC(xCov,gridTime,gridTS,Fs,pwrSpec,startLag,stopLag,5,'',PLOT_ACF,0.80,oneSideFq,dF);
            if PLOT_ACF
                if ranges(kl) < 2500 && ALTERNATE_COUNTER
                    title(['Hoover building lidar   ' datestr(gridTime(1),'yyyy-mm-dd HH:MM:SS') '  -  ' datestr(gridTime(end),'HH:MM:SS') '     Height ' num2str(round(ranges(kl))) ' m'])
                    print(['C:\Users\bcarroll\Documents\2022_DCFLUX\halo49\images\ACFplots\HooverBldgLidarACFplot_' datestr(gridTime(1),'yyyymmdd_HHMMSS') '_ht' num2str(round(ranges(kl)))],'-dpng')
                end
                if ALTERNATE_COUNTER==0
                    ALTERNATE_COUNTER=1;
                elseif ALTERNATE_COUNTER==1
                    ALTERNATE_COUNTER=0;
                end
            end
            crlbVarResult(kl) = nanmean(interp1(crlbSNRs,ideal_dCRLB,ib));
            varResult(kl) = nanvar(ts);
            nResult(kl) = length(ts);
            meanResult(kl) = mean(ts);
            skewResult(kl) = skewness(ts);
            kurtResult(kl) = kurtosis(ts);
            % Dissipation rate (e.g. Wulfmeyer et al 2024 Eq 18)
            if gridVarResult(kl).ACF_fitCoef < 0  % IF very little mixing, can end up with a negative ACF coef (negative slope in the lag fit range). Make those NaN, can't retrieve a reasonable dissipation rate.
                gridVarResult(kl).ACF_fitCoef = NaN;
            end
            epsilon(kl) = gridVarResult(kl).ACF_fitCoef^(2/3) / spdProf(kl);
            % Error in dissipation rate is from Wulfmeyer et al 2016 Appendix b.4, Eq A25.
            epsilonError(kl) = epsilon(kl) * sqrt(2.25*(gridVarResult(kl).ACF_fitCoef_confBounds/gridVarResult(kl).ACF_fitCoef)^2 + (WIND_SPD_UNC/spdProf(kl))^2);
            if CALC_EPSILON_AND_L_OC % O'Connor et al 2010 method
                res = calc_intEps_L(t,ts',spdProf(kl),[],[],[],[],[],dF);
                epsilon_OC(kl) = real(res.intEps);
                epsilonError_OC(kl) = real(res.epsE);
                integralScaleResult(kl) = res.L0v;
                integralScaleError(kl) = res.Le;
                integralScaleIsTime(kl) = res.integralScaleIsTime;
            end
            
        else  % ELSEIF fewer than minNumObs points in timeseries
            gridVarResult(kl) = nanVarResult;
        end
    end
    
    HaarResult = findCloudsAndBLH(snr2,decTime,10,ranges.*sind(nanmean(elev))-scannerHeight,3,params,params.BLH.HaarCloudHitThresh);
    gridVarResult = unpackStructArrays(gridVarResult); % swap to struct.val(index) format
    HaarResult = unpackStructArrays(HaarResult);
    
    % There are many ways to get variance that should be roughly equivalent. The
    % neatest approach for our current understanding/applications is to use only 
    % the ACF and lag0-lag1 for noise; no power spectrum.
    outVar = gridVarResult.totXcVar-gridVarResult.lag1Lag0Var; %-gridVarResult.acfNoiseVar; 
%     outVar = gridVarResult.totSpecVar-gridVarResult.acfNoiseVar;%gridVarResult.specNoiseVar;
%    outVar(gridVarResult.acfNoiseVar>1)=NaN;
    
    HFvar = HFvar - gridVarResult.acfNoiseVar.*pctElHFvar;
    HFvar(HFvar<0)=0;
    HFvar(gridVarResult.acfNoiseVar>1)=NaN;
    
    if displayGraphics
        lbl = sprintf('Top: Lidar wbSNR (dB) %s%s Bot: Lidar Vert Vel (m/s)',spc,spc);
        ttl = 'Lidar wbSNR w/ Cloud Base';
        
        figure(301)
        setCarboneColormap;
        clf
        %         ax(1) = subplot(211);
        ax(1) = subplot(2,4,[5,7]);
        % pcolor(decTime-jan0th,ranges,Vd');
        pcolor(decTime-jan0th,ranges,snr);
        shading flat;
        hold on
        set(gca,'Color',[0,0,0]);        % black backround
        %         datetick
        caxis(snrLims);
        cb  = colorbar;
        xlabel(cb,'dB');
        shading flat
        axis xy
        %
        ylim([0,yMax]);
        d = 0.5/60/24;
        xlim([nanmin(decTime-jan0th)-d,nanmax(decTime-jan0th)+d])
        set(gca,'Color',[0,0,0],'XTick',[nanmin(decTime-jan0th)-d:5/60/24:nanmax(decTime-jan0th)+d],'YTick',(0:yMax/numYTicks:yMax))
        datetickzoom('keeplimits','keepticks');
        
        ylabel('Height (m)');
        minTyear = datestr(startTime,'YYYY');
        %             lbl = sprintf('VOCALS RV Brown %s%s Initial day: %s (#%d)',spc,spc,datestr(useDate,'dd-mmm-yyyy'),julianDay);
        lbl = sprintf('%s %s%s%s Init day: %s (#%d)',titleStr,lilSpc,...
            'B: +Grd, G: -Grd, R: Var Yel:BLH Mag :CB ',lilSpc,datestr(testTime,'dd-mmm-yyyy'),julianDay);
        xlabel([{'Time (HH:MM) UTC'}]);
        ttl = 'Lidar wbSNR w/ Cloud Base';
        title([{ttl}]);
        ax(2) = subplot(2,4,[1,3]);
        pcolor(decTime-jan0th,ranges,vel);
        shading flat;
        hold on
        caxis([-5,5]);
        cb = colorbar;
        xlabel(cb,'ms^-^1');
        axis xy
        ylim([0,yMax]);
        d = 0.5/60/24;
        xlim([nanmin(decTime-jan0th)-d,nanmax(decTime-jan0th)+d]);
        set(gca,'Color',[0,0,0],'XTick',[nanmin(decTime-jan0th)-d:5/60/24:nanmax(decTime-jan0th)+d],'YTick',(0:yMax/numYTicks:yMax));
        datetickzoom('keeplimits','keepticks');
                
        ttl = 'Lidar W w/ Cloud Base';
        title([{ttl}]);
        ylabel('Height (m)');
        linkprop(ax(1:2),{'XLim'});
        
        
        a=1;
        varLims = [0,1];
        skewLims = [-2,2];
        snrLims = [-25,5];
        %         figure(2); clf
        ax(3) = subplot(2,4,4);
        
        %             plot(outVar,ranges,[0,1],blh*[1,1]); ylim([0,yMax]); xlim([0,1])
        %         ax(1)= gca;
        line(outVar,ranges,'Color','k');
        %             line([0,1],blh*[1,1],'Color','y');
%         for kl = 1: CI.numPeaks
%             if CI.numPtsInPeak(kl) / (CI.totPts-CI.numPtsNoCloud-CI.numBadPts) < 0.3 | CI.cloudFraction < 0.10 continue; end
%             line(varLims,CI.wtdPeakRange(kl)*[1,1],'color','m')
%             line(varLims,ranges(CI.minWidPos(kl))*[1,1],'color','m','linestyle','--')
%             line(varLims,ranges(CI.maxWidPos(kl))*[1,1],'color','m','linestyle','--');
%         end
        set(ax(3),'YTickLabel','');
        xlim(varLims); xlabel('Var (m^2/s^2)');
        ax(4) = axes('Position',get(ax(3),'Position'),...
            'YAxisLoc','right','XAxisLoc','top','Color','none',...
            'YColor','b','XColor','b');
        xlim(skewLims); xlabel('Skewness');
        ylabel('Height (m)');
        line(skewResult,ranges,'Color','b','Parent',ax(4));
        line([0,0],[0,max(ranges)],'Color','b','Parent',ax(4));
        
        ax(5) = subplot(2,4,8);
        %             plot(avgSnr,ranges,[-20,25],blh*[1,1],'y'); ylim([0,yMax]); xlim([-10,25]); hold on;
        plot(avgSnr,ranges,'b'); ylim([0,yMax]); xlim([-25,10]); hold on;
%         for kl = 1: CI.numPeaks
%             if CI.numPtsInPeak(kl) / (CI.totPts-CI.numPtsNoCloud-CI.numBadPts) < 0.3 | CI.cloudFraction < 0.10 continue; end
%             plot([-20,25],CI.wtdPeakRange(kl)*[1,1],'m',...
%                 [-20,25],ranges(CI.minWidPos(kl))*[1,1],'m--',...
%                 [-20,25],ranges(CI.maxWidPos(kl))*[1,1],'m--');
%         end
        set(ax(5),'YTickLabel','');
        xlim(snrLims); xlabel('SNR (dB)');
        % put ylabels on other side so does not interfere with colorbars
        ax(6) = axes('Position',get(ax(5),'Position'),...
            'YAxisLoc','right','XAxisLoc','top','Color','none',...
            'YColor','k','XColor','k');
        set(ax(6),'XTick',[]);
        ylabel('Height (m)');
        linkprop(ax,{'YLim'});
        ylim([0,yMax]);
        set(gcf,'color','white');
        suplabel(lbl,'x');
        [mov,oldMovDay] = outputGraphics(outFileStr,imageOutputDir,movieOutputDir,oneMovieFile,makeMovie,dumpGraphics,oldMovDay,mov,decTime);
        pause(0.1);
        a=1;
    end
    if outputData
        out.version = version;
        out.decTime = nanmean(decTime);
        out.lat = meanLat;
        out.lon = meanLon;
        out.avgSnr = avgSnr;
        out.varSnr = varSnr;
        out.avgBeta = avgBeta;
        out.wVar = outVar;
        out.wSkew = skewResult;
        out.wKurt = kurtResult;
        out.height = ranges + scannerHeight;
        out.n = nResult;
        out.fileNum = medianFileNum;
        out.HFvar=HFvar;
        out.rawVar = varResult;
        out.crlbVar = crlbVarResult;
        out.gridVarResult = gridVarResult;
        out.HaarResult=HaarResult;
        out.meanResult = meanResult;
        out.onlyUseCloudTopped = onlyUseCloudTopped;
%         if (onlyUseCloudTopped)
%             out.cloudBaseHeight = CI.meanPosNegGradRange - scannerHeight;
%             out.cloudBaseHeightVar = CI.varPosNedGradRange/CI.numPtsInPosNegGradPeak;
%             out.numBeamsInCloudBaseCalc = CI.numPtsInPosNegGradPeak;
%         else
%             out.cloudBaseHeight = CI.wtdPeakRange - scannerHeight;
%             out.cloudBaseHeightVar = CI.wtdPeakXSqrRange /CI.numPtsInPosNegGradPeak;
%             out.numBeamsInCloudBaseCalc = CI.numPtsInPeak;
%         end
        out.epsilon = epsilon;
        out.epsilonError = epsilonError;
        out.epsilon_OC = epsilon_OC;
        out.epsilonError_OC = epsilonError_OC;
        out.integralLengthScale = integralScaleResult;
        out.integralLengthScaleError = integralScaleError;
        out.integralScaleIsTime = integralScaleIsTime;
        outRecCount = outRecCount + 1;
        velStats(outRecCount) = out;
        outFilePending = 1;
    end
end

if outFilePending
%     outFileName = [outputDataDir,outFileStr,datestr(workDate,'_yyyymmdd'),'.mat'];
    outFileName = ['C:\Users\bcarroll\Documents\2022_DCFLUX\halo49\process output\profiles\velStats\',outFileStr,datestr(workDate,'_yyyymmdd'),'.mat'];
    fprintf('Writing to %s\n',outFileName);
    save(outFileName,'velStats');
end
 