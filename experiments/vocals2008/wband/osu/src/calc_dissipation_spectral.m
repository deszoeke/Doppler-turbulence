% Calculate cloud top TKE dissipation from Doppler velocity spectrum

%% set environment and constants
%cd /Users/sdeszoek/Data/cruises/VOCALS_2008/RHB/Scientific_analysis/programs/VOCALS2008_programs/wband
warning('off','MATLAB:interp1:NaNinY')
run('../read_parameters')
yyyy=sprintf('%04d',year);
% stday=318; % why missing time adjusted motion till 318? % 281;
% enday=336;
stday=318;
enday=336;
fs=3.5; %Hz
crosscomp=3/4;
kolmogorov=0.54;
factr=crosscomp/kolmogorov;

% w spectra window and FFT parameters used for dissipation (could pull out of loop)
%nwin=64;
%keep=2:8; % indices of spectra to keep - not good if detrended
%cut=10; % index at which to start regarding as noise
nwin=128;
keep=4:16; % indices of spectra to keep, start at 4 for detrended & windowed spectra, ?32 eats into noise floor?
cut=20; % index at which to start regarding as noise

F=(1:nwin/2)'/nwin*fs; % frequencies, Hz; 33: first is zero, last is Nyquist
dF=F(2)-F(1); % scalar
F53=F.^(5/3);
dt=10*60; % seconds per window (10 min})
di=floor(fs*dt); % samples per window

momentfile=dir([way_raw_data_wband '2008*MMCRMom.nc' ]);
h=ncread([way_raw_data_wband momentfile(1).name],'Heights',[1 1],[Inf 1]);

%% get motion status, read from Ken and Sergio's spreadsheet
motion_status
% 0 OK
% 1 no data
% 2 bias (<1 degree)
% 3 noisy
% 4 motion adj. off or failed

%% load 10-min ship-relative wind at cloud top interpolated from soundings
load('~/Data/cruises/VOCALS_2008/RHB/Scientific_analysis/programs/VOCALS2008_programs/wband/uv_cloudtop10.mat');
% ship maneuvers removed
% horizontal ship-relative speed
U=sqrt(utop10.^2+vtop10.^2); % m/s interpolated relative U from soundings and ship velocity

%% After best recalibration, Ken Moran says to subtract 1.72 dB from all recorded values.
% receiver gain is 1.72 higher than previously estimated, which reduces dBZ
% and dBmW signals. This is an increase in sensitivity. Note recalibration noise
% source error of +-0.6 dB.
dB_offset=-1.72; % dB
% define a noise level for the reflectivity based on Ken Moran's manual
radar_const=19.613; % dB_
analog_noise_signal=-120.74; % dBmW noise peak diagnosed by Simon
digital_noise_signal=-115.4;  % dBmW
noise_margin=+3.5; % dB
adhocthreshold=-43+dB_offset;
dBZnoiseavg=20*log10(h)+radar_const+analog_noise_signal;
dBZnoisefloor=20*log10(h)+radar_const+digital_noise_signal+noise_margin;
noisefloor=max(adhocthreshold,dBZnoisefloor);
%starter=find(base_time_yday>=310,1,'first');

%% loop length numbers
nhr=24*(enday-stday+1);
n10min=nhr*6;
%% allocate output data variables
count=zeros(n10min,1);
epsilon=NaN+zeros(n10min,1);
ydayw10=NaN+zeros(n10min,1);
sigv2=NaN+zeros(n10min,1);
spec=NaN+zeros(n10min,floor(nwin/2)+1);

%% loop over files by day, hour
ihr=0; % cumulative hour index
for iday=stday:enday;
    istatusday=find(iday==statusday); % index of motion status day (row)
    fprintf(1,'\n%3d ',iday)
    for idhr=0:23;
        istatushr=find(idhr==statushour); % index of motion status hour (col)
        fprintf(1,' %02d ',idhr)
        ihr=ihr+1;
        dddhh=sprintf('%03d%02d',iday,idhr);
        momentfile=dir([way_raw_data_wband '2008' dddhh '*MMCRMom.nc' ]);     
        kongfile=dir([way_raw_data_wband 'motion_adjT/' yyyy dddhh '*Kongsberg_adjT.txt']);
        if ~isempty(momentfile) && ~isempty(kongfile)
            countt=0; % number of data read so far
            [Z,vel]=deal(NaN+zeros(12600,120));
            %width=NaN+zeros(12600,120);
            for fi=1:size(momentfile,1) % handles multiple files
                %%prflname=[way_raw_data_wband momentfile(fi-1).name]; % previous file
                filename=[way_raw_data_wband momentfile(fi,:).name];
                
                % found no need to read the end of the previous file to get an integral minute
                % base_time in seconds since 1970-1-1 00:00:00
                base_time=ncread(filename,'base_time');
                % time_offset in seconds
                temp=ncread(filename,'time_offset');
                nt=length(temp);
                time_offset(countt+(1:nt),1)=temp;
                time(countt+(1:nt),1)=int32(time_offset(countt+(1:nt),1))+base_time;
                % base_time_mld in matlab datenumber
                base_time_mld=double(base_time)/86400 + datenum(1970,1,1,0,0,0);
                base_time_yday=base_time_mld-datenum(year,1,0,0,0,0);
%                 time_mld=double(time)/86400 + datenum(1970,1,1,0,0,0);
%                 time_yday=time_mld-datenum(year,1,0,0,0,0);
                time_yday_use(countt+(1:nt),1)=base_time_yday+time_offset(countt+(1:nt),1)/86400;
                
                % reflectivity matrix
                Z(countt+(1:nt),:)=ncread(filename,'Reflectivity')'+dB_offset; % new matlab2011a ncread
                %+dB_offset 2013-04-28 SPdeS (1.72 dB more liberal cloud detection was used for the 2011 paper);
                vel(countt+(1:nt),:)=ncread(filename,'MeanDopplerVelocity')';
                %width(countt+(1:nt),:)=ncread(filename,'SpectralWidth')';
                countt=countt+nt; % number of data read so far
            end % multiple files
            
            % Read Kongsberg motion compensation file
            countk=0;
            for ki=1:length(kongfile)
                kongfilename=[way_raw_data_wband 'motion_adjT/' kongfile(1).name]; % _adjT data starts on day 318
                [ktimetmp,kwtmp]=read_kongsberg(kongfilename);
                nk=length(ktimetmp);
                kongtime(countk+(1:nk),1)=ktimetmp;
                kongw(countk+(1:nk),1)=kwtmp;
                countk=countk+nk;
            end
            
            % synchronize Kongsberg and radar time (yearday)
            [junk,ii,jj]=unique(kongtime);
            kongw4radar=interp1(kongtime(ii),kongw(ii),time_yday_use(1:countt));
                      
            % cutoff below noise floor
            mask=Z>=repmat(noisefloor',[length(Z) 1]);
            Z(~mask)=NaN;
            vel(~mask)=NaN;
            
            % ship-motion-compensated vertical velocity +towards
            w=vel(1:countt,:)+repmat(kongw4radar(1:countt),[1,length(h)]);

            %width(~mask)=NaN;
            % cutoff below 170 m range gate (clutter, near field contamination)
            Z(:,h<170)=NaN;
            w(:,h<170)=NaN;
            %width(:,h<170)=NaN;
            Z=Z(1:nt,:); % truncate
            w=w(1:nt,:);
            %width=width(1:nt,:);
            
            % calculate cloud-top dissipation from inertial range of spectra            
            % procedure from diag_dissipation.m for calculating dissipation at cloud top...
            
            % loop through, interp and calc spectra for 10-min chunks
            for ind=1:round(nt/di) % index of the 10 minute interval in the hour
                bigind=24*6*(iday-stday)+6*idhr+ind;
                i10=(ind-1)*di+1; % first index in the 10 minute interval
                ydayw10(bigind,1)=iday+(idhr+(ind-1)/6)/24; % current time of this 10-minute window
                % Utop=nanmean(U(yday10>=ydayw10(bigind) & yday10<ydayw10(bigind)+10/60/24));
                Utop=double(interp1(yday10,U,ydayw10(bigind),'nearest')); % interp ship-rel wind speed to the current time
                % what is the max coverage of nonnoise returns? --> nmax
                ii10=i10:min(i10+di-1,nt);
                if length(ii10)<di/2; continue; end % break from loop if too few obs.
                nw=sum(isfinite(w(ii10,:)));
                [nmax,knmax]=max(nw);
                % choose the highest level that has >0.5*nmax of this coverage --> ik
                ik=find(nw>0.5*nmax,1,'last');
                if ~isempty(ik) && ik<100
                    w10=w(ii10,ik-1:ik+1);
                    iwt=isfinite(w10);
                    % take highest finite w within a level of this 50%ile cloud top
                    wtop=w10(:,1); % w is positive-down
                    wtop(iwt(:,2))=w10(iwt(:,2),2); % overwrites
                    wtop(iwt(:,3))=w10(iwt(:,3),3); % overwrites if there's a higher one
                    % filter out extreme values
                    wtop(abs(wtop-nanmean(wtop))>3*nanstd(wtop))=NaN;
                    % interpolate missing values
                    isn=isnan(wtop);
                    wtop(isn)=interp1(time_yday_use(ii10(~isn)),wtop(~isn),time_yday_use(isn));
                    count(bigind,1)=sum(isfinite(wtop));
                    if count(bigind,1)>di/2;
                        %[S,F]=pwelch(wtop(isfinite(wtop)),nwin,floor(nwin/2),nwin,fs); % uses hamming window 0.54-0.46cos % 0.036 s elapsed
                        %[Sm,Fm]=pwelch(wtop(isfinite(wtop)),hann(nwin),floor(nwin/2),nwin,fs); % hanning cos window (makes no difference)
                        [S,F]=fast_psd(wtop(isfinite(wtop)),nwin,fs); % 2x faster, Hanning window and detrends before FFT
                        % low power in Nyquist frequency (end) could be from sampling out of phase with cycles
                        % probably has to do with discrete time series unable to sample phase and quadrature phase
                        % Ignore Nyquist-frequency power estimate!
                        %spec(bigind,:)=S';
                        Snoise=min(S(cut:end-1));
                        vls=factr.*(2*pi/Utop)^(2/3).*mean(F53(keep).*(S(keep)-Snoise)); % dissipation ^2/3
                        epsilon(bigind,1)=vls.^1.5; % dissipation
                        sigv2(bigind,1)=dF*sum(S(2:end)-Snoise); % calculate variance from the resolved nonnoise spectrum
%                         % Emily Shroyer's method of fitting
%                         % f^-5/3 and a noise floor 
%                         % does not seem to work here
%                         Form=[ones(size(F(4:end-1))),F(4:end-1).^(-5/3)];
%                         FIT2=lsqlin(Form,S(4:end-1),[],[],[],[],[0 0],[  ],[epsilon(bigind,1) Snoise]);
%                         eeps=2*pi/Utop*(factr*FIT2(2))^(3/2);
%                         % test plot spectrum and fits
%                         loglog(F,S,F,1./F53*4/3*0.54*(Utop/2/pi*epsilon(bigind,1))^(2/3)+Snoise,...
%                             F,1./F53*4/3*0.54*(Utop/2/pi*eeps)^(2/3)+FIT2(1))
%                         hold on
%                         loglog(F,FIT2(2)./F53+FIT2(1),'m')
                    end % there are spectra
                end % there is 10 min data
            end % 10 min
        end % there is hourly data
    end % hr
end % day

% truncate data
ydayw10=ydayw10(1:bigind);
sigv2=sigv2(1:bigind);
epsilon=epsilon(1:bigind);
count=count(1:bigind);
spec=spec(1:bigind,:);

%% QC for motion status
% motion_status; % read above
% 0 OK
% 1 no data
% 2 bias (<1 degree)
% 3 noisy
% 4 motion adj. off or failed
st=Status';
statusyday=repmat(statusday,[1 24])+repmat(statushour/24,[23 1]);
statusydayt=statusyday';
% plot(statusydayt(:),st(:))
% convert hourly to 10-min
a=repmat(st(:),[1 6])';
status10=a(:);
a=(repmat(statusydayt(:),[1 6])+repmat((0:5)*10/60/24,[length(statusydayt(:)) 1]))';
statusyday10=a(:);
statusint=interp1(statusyday10,status10,ydayw10+1/60/24,'nearest');
% ADD status 5 to statusint: no dissipation retrieval
statusint(isnan(statusint))=5;

%% QC for enough points in the spectrum
isenuf=count>1800;

%% QC for missing motion data and motion correction failure, 85% data available
%  and, a posteriori, variance <0.4 m^2/s^2
isgdtmp=isenuf & (statusint==0 | statusint==2 | statusint==3); % bias and noisy motion correction seem to be OK for spectra
% nanstd(sigv2(isgdtmp))*3 ~ 0.4
isgd=isgdtmp & sigv2<0.4;
sp=spec(isgd,:);
ep=epsilon(isgd);
yd=ydayw10(isgd);


%% save data
save cloudtopturb.mat ydayw10 F spec epsilon sigv2 count statusint

%% recover noise for each spectrum
snoise=min(spec(:,cut:end-1),[],2);

%%% plot

%% test QC procedure
plot(ydayw10,statusint/10,'y.')
hold on
plot(ydayw10(isgd),sigv2(isgd),'r.')
plot(ydayw10,sigv2,'b')
plot(ydayw10,epsilon/1e-2,'k')
plot(ydayw10(isgd),epsilon(isgd)/1e-2,'m.')
set(gca,'ylim',[0 1])

%% scatter plot epsilon and power
plot(epsilon(isgd),sigv2(isgd),'.')
hold on
set(gca,'yscale','log','xscale','log')
axis([0 0.01 0 1])
plot([1.7e-3 1.7e-2].^1.5,[2e-2 2e-1],'r') % epsilon=power^1.5 dependence
xlabel('dissipation')
ylabel('<w''w''>')

%% time series
subplot(2,1,1)
x=ydayw10(isfinite(ydayw10)); y=epsilon(isfinite(ydayw10));
isx=isgdtmp(isfinite(ydayw10)) & sigv2(isfinite(ydayw10))<0.4;
% y(~isx)=NaN;
% plot(x-5/24,y)
y(~isx)=-.01e-3;
area(x-5/24,y,'facecolor',[0.5 0.6 0.7],'edgecolor','none')
hold on
plot(ydayw10(isgd)-5/24,epsilon(isgd),'.','markersize',4)
plot(ydayw10-5/24,5e-4*max(0,-cos(2*pi*(ydayw10-5/24))))
x6=x(1:6:end);
y6=all(~isx([1:6:end; 2:6:end; 3:6:end; 4:6:end; 5:6:end; 6:6:end]))'
plot(x6-5/24,.7e-3*y6-.1e-3,'r.')
set(gca,'tickdir','out','ylim',[0 3e-3],'fontsize',14,'xminortick','on')
xlabel('2008 yearday')
ylabel('TKE dissipation m^2 s^{-3}')
print -depsc CTdissipation/dissip_timeseries.eps

%% QC outcome displayed diurnally
plot(mod(ydayw10*24-5,24),epsilon,'rx') % all points
hold on
plot(mod(ydayw10(isgd)*24-5,24),epsilon(isgd),'.') % only good points
% histogram
subplot(2,1,1)
hist(log10(epsilon),-5.5:0.05:-0.5)
title('all')
set(gca,'xlim',[-5.5 -0.5])
subplot(2,1,2)
hist(log10(ep),-5.5:0.05:-0.5)
set(gca,'xlim',[-5.5 -0.5])
title('QCd')
xlabel('log_{10}(\epsilon)')
% dissipation is lognormally distributed
% QCing based on count and motion correction has little effect on its distribution

%% shade spectra and plot 10x<w'w'>
ii=isfinite(sigv2);
pcolor(ydayw10(ii),F,log10(spec(ii,:)')); shading flat
hold on
plot(ydayw10,10*sigv2,'k')

%% sort the spectra by power in low, high freq tails
clf
subplot(2,1,1)
[ssp,iord]=sort(sp(:,2));
pcolor(log10(sp(iord,:))'); shading flat
hold on
plot(1e4*ep(iord),'k')
set(gca,'xlim',[1 length(iord)],'ylim',[0 33])

subplot(2,1,2)
[ssp,iord]=sort(sp(:,end-1));
pcolor(log10(sp(iord,:))'); shading flat
hold on
plot(1e4*ep(iord),'k')
set(gca,'xlim',[1 length(iord)],'ylim',[0 33])

%% divide spectra into modes with an EOF analysis
ii=isfinite(sigv2);
X=log10(spec(ii,2:end-1))-repmat(mean(log10(spec(ii,2:end-1))),[sum(ii) 1]); % de-mean log10 power series
[u,s,v]=svd(X);

subplot(2,2,1)
pcolor(X')'; shading flat
title('full matrix')
subplot(2,2,2)
p=zeros(size(s)); p(1,1)=s(1,1);
pcolor((u*p*v')'); shading flat
title('mode 1')
hold on
plot(v(:,1)/max(v(:,1))*1e3,1:31,'color','k','linewidth',1.4)
subplot(2,2,3)
p=zeros(size(s)); p(2,2)=s(2,2);
pcolor((u*p*v')'); shading flat
title('mode 2')
hold on
plot(v(:,2)/max(abs(v(:,2)))*5e2+5e2,1:31,'color','k','linewidth',1.4)
subplot(2,2,4)
p=zeros(size(s)); p(3,3)=s(3,3);
pcolor((u*p*v')'); shading flat
title('mode 3')
hold on
plot(v(:,3)/max(abs(v(:,3)))*5e2+5e2,1:31,'color','k','linewidth',1.4)

%% composite diurnal cycle of dissipation on local hour: mean, median histogram
indx=bin2(0:24,mod(yd*24-5,24)); % local hour index 1:24; 1 corresp to 0-1 local
ledges=-5.5:0.1:-0.5;
edges=0:1e-4:1e-2;
for ih=1:24;
    ii=indx==ih;
    epse(ih)=nanmean(ep(ii));
    epsm(ih,:)=quantile(ep(ii),[.25 .50 .75]);
    % histogram
    hlep(ih,:)=histc(log10(ep(ii)),ledges);
    hep(ih,:)=histc(ep(ii),edges);
    % histogram with median diurnal cycle removed
    y=ep(ii)-epsm(ih,2);
    hdlep(ih,:)=histc(log10(y),ledges);
    hdep(ih,:)=histc(y,edges-2e-3);
end

clf
plot(mod(yd*24-5,24),ep*1e3,'.','color',[0.7 0.7 0.7])
hold on
plot(0.5:23.5,epse*1e3,'ko') % local hour
lh=plot(-0.5:24.5,epsm([end 1:end 1],:)*1e3,'k'); % local hour
set(lh(2),'linewidth',2)
set(gca,'fontsize',16,'xlim',[0 24],'xtick',0:6:24,'xticklabel',[0:6:18 0])
ylabel('TKE dissipation (10^{-3} m^2 s^{-3})')
xlabel('local hour')
set(gca,'ylim',[0 1.5])
print -depsc CTdissipation/dissipation_diel.eps

colormap(1-gray(16))
% log scale diurnal histogram
subplot(2,1,1)
pcolor(0:36,ledges,hlep([1:end 1:end/2+1],:)'); shading flat
colorbar
set(gca,'ylim',[-5 -2],'xtick',[0:6:36],'xticklabel',[0:6:18 0:6:12],'tickdir','out')
subplot(2,1,2)
pcolor(0:36,ledges,hdlep([1:end 1:end/2+1],:)'); shading flat
colorbar
set(gca,'ylim',[-5 -2],'xtick',[0:6:36],'xticklabel',[0:6:18 0:6:12],'tickdir','out')

% linear scale diurnal histogram
subplot(2,1,1)
pcolor(0:36,edges,log2(hep([1:end 1:end/2+1],:))'); shading flat
colorbar
hold on
plot(-0.5:36.5,epsm([end 1:end 1:end/2+1],2),'g')
set(gca,'ylim',[0 2e-3],'xtick',[0:6:36],'xticklabel',[0:6:18 0:6:12],'tickdir','out')
subplot(2,1,2)
pcolor(0:36,edges-2e-3,hdep([1:end 1:end/2+1],:)'); shading flat
colorbar
set(gca,'ylim',[-0.5e-3 2e-3],'xtick',[0:6:36],'xticklabel',[0:6:18 0:6:12],'tickdir','out')

%% log2-frequency diurnal histogram
clf
pcolor(0:36,edges*1e3,log2(hep([1:end 1:end/2+1],:))'); shading flat
a=b2rcolormap(20);
caxis([-0.5 5])
cm=colormap([[1 1 1]; a([11 11:19],:)]);
hc=colorbar;
set(hc,'ylim',[-0.1 4.6],'ytick',[0 1:0.5:4.5],'yticklabel',[1 2 3 4 6 8 12 16 23] ,'fontsize',14)
hold on
plot(mod(yd*24-5,24),ep*1e3,'k.','markersize',3)
plot(24+mod(yd*24-5,24),ep*1e3,'k.','markersize',3)
hml=plot(-0.5:36.5,epsm([end 1:end 1:end/2+1],:)*1e3,'color',[0 0 0]);
set(hml(2),'linewidth',2)
plot(0.5:35.5,epse([1:end 1:end/2])*1e3,'ko') % local hour
set(gca,'ylim',[0 2],'xlim',[0 36],'xtick',[0:6:36],'xticklabel',[0:6:18 0:6:12],...
    'tickdir','out','fontsize',16)
xlabel('local hour')
ylabel('cloud top dissipation (10^{-3} m^2 s^{-3})')
print -depsc CTdissipation/dissipation_diel_shade.eps

% de-meaning by shifting the peak of the distribution takes out recognizable
% cycle of the peak but broadens the distribution of the tails of the distribution of epsilon.
% In the raw series it may be bounded below by noise.

%% distribution shifts over the diurnal cycle
% start at 17 local, 24
b2rcolormap(25);
subplot(2,1,1)
area(ledges,hlep([18:end 1:17],:)')
set(gca,'xlim',[-5 -2],'fontsize',16)
xlabel('log_{10}\epsilon')
subplot(2,1,2)
area(edges,hep([18:end 1:17],:)')
set(gca,'xlim',[0 5e-3],'fontsize',16)
xlabel('\epsilon')

% shifed to remove cycle of diurnal median
subplot(2,1,1)
area(ledges,hdlep([18:end 1:17],:)')
set(gca,'xlim',[-6 -2],'fontsize',16)
xlabel('log_{10}\epsilon')
title('diurnal median removed')
subplot(2,1,2)
area(edges,hdep([18:end 1:17],:)')
set(gca,'xlim',[1e-3 5e-3],'fontsize',16)
xlabel('\epsilon')

