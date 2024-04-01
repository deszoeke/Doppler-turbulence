% DEW KM read in W-Band nc moment data and plot
% range gates are 25 m
% timing is every 0.25 sec
%Version_3 reads kongsberg motion data and corrects vertical velocity
%   cd c:\data\cwf\matlabstf\cwf\VOCALS2008_programs_cwf\wband
%The read_parameters_VOCALSz.m file sets the paths for the radar data
% 96  W-band radar files
% 97  way_raw_data_wband=[big_dir '\RHB\radar\wband\Raw\'];
% 98  way_raw_images_wband=[big_dir '\RHB\radar\wband\Raw_images\'];
% You must go to lines 97-98 and set the path for your data archive
%close all
clear
read_parameters;
newfigures=1;
motionread=1;%set to zero to skip motion correction

year=2008;
% ddds = input('Input yearday  ','s');
% shrs = input('Input start hr  ','s');
% ehrs = input('Input end hr  ','s');
ddds='325'; shrs='12'; ehrs='12';
ddd = str2num(ddds);
shr = str2num(shrs);
ehr = str2num(ehrs);

%eval(['cd ', way_raw_data_wband]);%eval(['cd ', 'C:\data\cwf\data\vocals\']);
save dew1 ddds ddd shrs ehrs shr ehr
%ncload('C:\Data\VOCALS\Data\wband\NetCDFRadarmfiles\WbandMom.nc')  % First rain shower
%20082821346MMCRMom.nc  20082821400MMCRMom.nc% 2nd rain shower

ht = (1:1:120)*.025; % 1st guess range gates in km
fn = [way_raw_data_wband '20082821400MMCRMom.nc'];
fn(41:43) = ddds;

cldtop=[];
for jam = shr:ehr
    save dew2 jam
    %close all
    clear hh y z dt jd
    read_parameters;
    load dew1
    load dew2

    %ncload('C:\Data\VOCALS\Data\wband\NetCDFRadarmfiles\WbandMom.nc')  % First rain shower
    %20082821346MMCRMom.nc  20082821400MMCRMom.nc% 2nd rain shower
    
    hr=sprintf('%02d',jam);
    fn = [way_raw_data_wband '2008' ddds hr '00MMCRMom.nc'];
    if exist(fn,'file') % Test to see if file exists
        ncload(fn)
    else
        disp(['File not found: ' fn])
        continue % to next hour
    end
    
    tt = time_offset + base_time;
    % calculate time from time_offset and base_time
    y=tt/(60*60*24);
    z=y+datenum([1970,1,1,0,0,0]);
    dt = datestr(z, 31);    %2008-10-07 16:35:53
    [yr,month,day,hour,minutes,secs]=datevec(z);
    jd=datenum(yr,month,day,hour,minutes,secs)-datenum(yr-1,12,31);
    hh =(jd-(fix(jd(1))))*24; % subtract off JD and convert to an hour
    
    ht = Heights(1,:)/1e3; % range gates in km
    Ref = Reflectivity'; % Ref(height, time)
    SW = SpectralWidth';
    DVa = MeanDopplerVelocity';
    SNR = SignalToNoiseRatio';
    DV=NaN+DVa;
    
    if motionread %correct doppler for
        %         read_Wband_motion                   % wasn't working
        %         DV_mot=interp1(hhk,vkons(9,:),hh);
        %         ii=find(isnan(DV_mot));DV_mot(ii)=0;
        %         for im=1:length(ht);
        %             DV(im,:)=DVa(im,:)+DV_mot;
        %         end;
        % Read Kongsberg motion compensation file
        kongfile=dir([way_raw_data_wband 'motion_adjT/' sprintf('%04i%03i%02i',[year,ddd,jam]) '*Kongsberg_adjT.txt']);
        if isempty(kongfile)
            kongflag=1; % error no file
            DV=DVa;
        elseif length(kongfile)>1
            kongflag=2; % error multiple files
            DV=DVa;
        else
            kongflag=0;
            kongfilename=[way_raw_data_wband 'motion_adjT/' kongfile(1).name];
            [kongtime,kongw,kongpitch,kongroll]=read_kongsberg(kongfilename);
            
            % synchronize Kongsberg and radar time (yearday)
            base_time_mld=double(base_time)/86400 + datenum(1970,1,1,0,0,0);
            base_time_yday=base_time_mld-datenum(year,1,0,0,0,0);
            [junk,ii,jj]=unique(kongtime);
            DV_mot=interp1(kongtime(ii),kongw(ii),base_time_yday+time_offset/86400)';
            DV=DVa+repmat(DV_mot,[120 1]);
        end
    else
        DV=DVa;
    end;
    
    clear Reflectivity SpectralWidth MeanDopplerVelocity SignalToNoiseRatio
    for k=1:length(Ref);
        kk=find(Ref(:,k)+39-20*log10(ht')>8 & ht'>.5 & ht'<3 & Ref(:,k)>-30);
        if ~isempty(kk)
            htt(k)=ht(kk(length(kk)));
            reft(k)=mean(10.^(Ref(kk,k)/10));
            reft2(k)=mean(10.^((Ref(kk,k)+39-20*log10(ht(kk)'))/10));
            wvel(k)=mean(DV(kk,k));
            nnt(k)=length(kk);
        else
            htt(k)=NaN;
            reft(k)=NaN;
            reft2(k)=NaN;
            wvel(k)=NaN;
            nnt(k)=0;
        end;
    end;
    y = prctile(htt,[10 15 50 85 90]);
    kk=find(htt>=y(1) & htt<=y(5));
    httb=nanmean(htt(kk));
    reftb=10*log10(nanmean(reft(kk)));
    reft2b=10*log10(nanmean(reft2(kk)));
    wvelb=nanmean(wvel(kk));
    nntb=sum(nnt);
    sigh=sqrt(nanmean(htt(kk).^2)-httb.^2);
    disp(['Hour= ' num2str((hh(1))) '  Mean cloud top ht= ' num2str(httb) '  std cloud top height = ' num2str(sigh)])
    disp(['Mean dBZ= ' num2str(reftb) '   Mean dBZ-noisedBZ= ' num2str(reft2b) '   Mean Vel= ' num2str(wvelb)])
    cldtop=[cldtop' [hh(1) httb sigh (y(4)-y(3))/2 reftb reft2b wvelb nntb]']';
    
    if newfigures; figure(1); dock; end; clf
    axlim=[12 24 0 2];
    subplot(3,1,1),imagesc((hh-jam)*60,ht, Ref(1:120,:),[-60, 20]);
    axis xy; axis(axlim)
    ylabel('Ht (km)')
    title({[cruise yyyy ' W-Band'' YD ' ddds ' Hr ' hr]})
    h = colorbar;
    axes(h); ylabel('Ref dBZ')
%     subplot(4,1,2),imagesc((hh-jam)*60,ht, SNR(1:120,:),[-25, 20]);  %
%     axis xy
%     ylabel('Ht (km)')
%     %axis([15 15.25 0 3])
%     h = colorbar;
%     axes(h)
%     ylabel('SNR dB')
    subplot(3,1,2),imagesc((hh-jam)*60,ht, DV(1:120,:),[-5 5]);  %
    axis xy; axis(axlim)
    ylabel('Ht (km)')
    h = colorbar;
    axes(h); ylabel('DopVel (m/s)')   
    subplot(3,1,3),imagesc((hh-jam)*60,ht, SW(1:120,:),[0, 5]);  %
    axis xy; axis(axlim)
    ylabel('Ht (km)')
    xlabel(sprintf('Minute (UTC hour %02d)',jam))
    h = colorbar;
    axes(h); ylabel('SpecWidth (m/s)')
    %print('-dpng','Wband_Ref_DV_SW325_12_12-24.png');     
%     subplot(4,1,4),imagesc((hh-jam)*60,ht, 10*log10(SW./abs(DV)),[-10 20]);  %
%     axis xy; axis(axlim)
%     ylabel('Ht (km)')
%     h = colorbar;
%     axes(h); ylabel('DoppVel Width/Mean (dB)')

end % hr loop

if false % SPdeS 2009/09/17 diagnostics of strange spectral width 2008 325 12:15
    % Fairall sent spectral velocity pics from a different time that shows the radar
    % sometimes has a velocity side lobe problem.
    s=find((hh-jam)*60>15,1,'first'); % minute 15
    % subplot(1,2,1)
    plot(SW(:,s+[0:20*3.5-1]),repmat(ht',[1 20*3.5])+.025*rand(size(ht,2),20*3.5)-.0125,'k.','markersize',4)
    axis([0 4 0 2.5])
    title('VOCALS 325 12:15:00-12:15:20')
    xlabel('Spectral width (m/s)')
    ylabel('height (km)')
    %print -dpng Wband_ht_width_325_12_15_00-20.png
    % s=find((hh-jam)*60>18,1,'first'); % minute 18
    % subplot(1,2,2)
    % plot(SW(:,s+[0:20*3.5-1]),repmat(ht',[1 20*3.5])+.025*rand(size(ht,2),20*3.5)-.0125,'.','markersize',4)
    % axis([0 4 0 2.5])
    % title('VOCALS 325 12:18:00-12:18:20')
    % xlabel('Spectral width (m/s)')
    
    z=ht>=1.1 & ht<=1.4;
    plot(SW(z,s+[0:60*3.5-1]),DV(z,s+[0:60*3.5-1]),'b.','markersize',4)
    hold on
    plot(SW(z,s+[0:60*3.5-1])+DV(z,s+[0:60*3.5-1]),DV(z,s+[0:60*3.5-1]),'k.','markersize',4)
    % SW vs DV reveals nothing else.
end

%%%%%%%%%%%%%%  Some diagnostics CWF added on the vocals cruise
ceilo=1;
if ceilo
    zz=load([way_proc_data_ceilo, cruise, yyyy, 'ceilo_h_a.txt']);
    jdh=zz(:,1)+.5/24;%1    Julian Day (2003), time at the beginning of 10-min ave
    allda=zz(:,2);%number of obs in 10-min
    nocld=zz(:,3);%no clouds
    onecld=zz(:,4);%one cloud
    mltcld=zz(:,5);%multiple clouds
    obscure=zz(:,6);%obscured
    pobs=zz(:,7);%partially obscured
    clrf=zz(:,8);%clear fraction
    ceilof=zz(:,9);%cloudy fraction
    ceilofobs=zz(:,10);%cloudy fraction, including obscured
    cloudhgt=zz(:,11);%median cloud base height
    hgts1=zz(:,12);%15% clouds lower height
    hgts2=zz(:,13);%85% clouds lower height
end

if length(cldtop)>22;
    %ddd=327;ddds=num2str(ddd);cldtop=load(['cldtop_' ddds '.txt']);
    
    jdt=ddd+cldtop(:,1)/24+.5/24;
    ii=find(cldtop(:,8)<500);cldtop(ii,2:7)=NaN;
    
    figure;
    subplot(2,1,1);plot(jdt,cldtop(:,2),'-o',jdt,cldtop(:,2)+cldtop(:,3),'--',jdt,cldtop(:,2)-cldtop(:,3),'--',jdt,cldtop(:,2)-cldtop(:,8)/12600*25/1e3,'x-');
    if ceilo;hold on;plot(jdh,cloudhgt/1000,'-d',jdh,hgts1/1000,'.',jdh,hgts2/1000,'.');axis([ddd ddd+1 .5 1.8]);end;
    Ylabel('Cloud top/base height (km)');
    subplot(2,1,2);plot(jdt,cldtop(:,5),'-o',jdt,10*cldtop(:,7),'-x',[ddd ddd+1],[0 0],'-' ,[ddd ddd+1],[-17 -17],'--' );
    Ylabel('dBZ and 10*W (m/s)');
    xlabel('Julian Day (2008)');
    ii=find(floor(jdh)==ddd);
    figure;plot(jdt,cldtop(:,2)-cloudhgt(ii)/1000);axis([ddd ddd+1 0 .5]);
    Ylabel('Cloud Thickness (m)');
    xlabel('Julian Day (2008)');
    
    save ([ddds '.txt'],'cldtop','-tabs','-ascii')
end;

[nrange ntime]=size(Ref);
H=repmat(ht',[1 ntime]); % spdes
% detect=Ref+39-20*log10(H)>8 & H>.1 & H<3; % detect threshold 8 dBZ above noise
detect=Ref+39-20*log10(H)>4 & H>0.250 & H<1.950; % detect threshold 4 dBZ above noise

% exclude low-Z returns using detect threshold
xRef=Ref; % dBZ
xDV=DV;
xSW=SW;
xZZ=10.^(Ref/10); % Z (mm6/m3 radar units = 10^-18 m3 SI units)
xRef(~detect)=NaN;
xDV(~detect)=NaN;
xSW(~detect)=NaN;
xZZ(~detect)=NaN;

xRefm=nanmean(xRef,2); % average only returns when drizzle is present
xZZm=nanmean(xZZ,2);
xDVm=nanmean(xDV,2);
xSWm=nanmean(xSW,2);
detectm=mean(detect,2);
% To convert unconditional time average into a drizzle-only average divide by drizzle_fraction=detectm.

mm=sum(detect,2)<5;
xRefm(mm)=NaN; %  bug fix, was xRevm
xDVm(mm)=NaN;
xZZm(mm)=NaN;
%xSWm(mm)=NaN; % do not filter xSWm more, bc. want to get a floor on the width above the cloud

dBZbar=10*log10(xZZm);
dBZbar(dBZbar<-60)=NaN;

if newfigures; figure(2); dock; end; clf
subplot(1,2,1);
plot(xRefm,ht,'color',[0 0 1])                            % detect>5 mean
hold on
plot(mean(Ref,2),ht,'color',[0 0.5 0])                    % unconditional mean
plot(dBZbar,ht,'--','color',[1 0 0])                      % detect>5, dBZ>-60 mean (avg linear Z)
plot(max(-60,-39+20*log10(ht')),ht,'color',[0 0.75 0.75]) % noise floor
plot([-17 -17],[0 3],'--','color',[0.75 0 0.75])          % threshold
ylabel('Height (m)'); xlabel('dBZ'); %axis([-60 20 0 2]);
title 'reflectivity'
subplot(1,2,2);
plot(detectm,ht);
xlabel('Return Fraction');%axis([0 1 0 2]);

% Retrieve the lognormal drizzle parameters, assuming the detected velocities
% are fall velocities.
[SigX, ModalRadius, N]=drizzle_lognorm_param(xDV,xSW,xZZ*1e-18); % Z[mm6/mm3]-->[m3]
drizzledetect=xRef>-5 & xDV>0.3 & xDV<3.0;
% dr=0.01e-4;
% edges=(0:dr:2e-4)';
logr=(-11:.01:-8.5)';
edges=exp(logr);
df=diff(edges(1:end)); df(end+1)=df(end);
n_cloud=histc(ModalRadius(~drizzledetect),edges)./df;
n_drizzle=histc(ModalRadius(drizzledetect),edges)./df;
% integral-scale type average of modal radius
modal_radius_cloud=dr*sum(n_cloud)/n_cloud(1);
modal_radius_drizzle=dr*sum(n_drizzle)/n_drizzle(1); % cloud vs. drizzle condition makes negligible difference
clf
loglog(edges,n_cloud./sum(n_cloud),'.')
hold on
loglog(edges,n_drizzle./sum(n_drizzle),'rx')
legend('cloud','drizzle')
xlabel('modal radius (m)')
ylabel('normalized frequency')
%print -dill lognormal_modalradius.ai % subsequently edited in AIllustrator
ii=edges>45e-6 & n_drizzle>0;
pf=polyfit(log(edges(ii)),log(n_drizzle(ii)./df(ii)),2);
sigxest=sqrt(-1/pf(1)); % estimate of sigma(log(r))
xmax=-pf(2)/(2*pf(1)); % corresp to 56.4 microns


% average above cloud
ha=ht>2; % height km above which to average
wmm=mean(mean(DV(ha,:),2)); % avg in time first -  Doppler velocity
wdm=mean(mean(SW(ha,:),2)); %                      Spectral velocity width

% noise-filtered
DVbar=xDVm-wmm;
dBZ=xRefm;
Widthbar=xSWm;
wdbm=nanmean(Widthbar(ha,:));

% noisier, average weak Z returns too
DVbar2=mean(DV,2)-wmm;
dBZ2=mean(Ref,2);
ZZ=mean(10.^(Ref/10),2);
dBZbar2=10*log10(ZZ);
Widthbar2=mean(SW,2);
wdbm2=nanmean(Widthbar2(ha,:)); % spectral width of clear air (noise SW)

if newfigures; figure(3); dock; end; clf
subplot(1,2,1);
plot(DVbar,ht,'color',[0 0 1])
hold on
plot(DVbar2,ht,'--','color',[0 0.5 0])
plot([0 0],[0 3],'color',[1 0 0]);
ylabel('Height (m)');xlabel('Mean Velocity (m/s)');%axis([-1 5 0 2]);
subplot(1,2,2);
plot(Widthbar,ht,'color',[0 0 1])
hold on
plot(Widthbar2,ht,'--','color',[0 0.5 0])
plot([0 0],[0 3],'color',[1 0 0]);
xlabel('Doppler Width (m/s)'); axis([0 1.5 0 3]);

%Drizzle microphysics parameterization from Frisch et al 1995
A=1.2e-4; % s  - Gossard et al 1990 constants for linear drop fall regime (drizzle)
B=1e-5;   % m

%sigv^2 is (spectral width)^2 - (noise spectral width)^2 ; this is not eqn.12!
sigvsq=Widthbar.^2-wdbm^2;
%sigx=sqrt(Widthbar.^2-wdbm^2)./DVbar; % approx...
sigx=sqrt(log(1+sigvsq./(DVbar+B/A).^2));
sigx(sigx>1)=NaN; DVbar(sigx>1)=NaN;
yz=dBZ-132+26*sigx.^2;
Ql=(A*DVbar).^(-3).*10.^(yz/10);
Ql(Ql>10)=NaN;%g/m^3
Fq=Ql.*((DVbar+B/A).*exp(-3*sigx.^2)-B/A)*1e-6*1000*3600;%mm/hr
r0=1e3*(A*DVbar+B).*exp(-6.5*sigx.^2);%drizzle radius in mm

% sigx2=sqrt(Widthbar2.^2-wdm^2)./DVbar2;
sigx2=sqrt(log(1+(Widthbar2.^2-wdm^2)./(DVbar2+B/A).^2));
sigx2(sigx2>1)=NaN; DVbar2(sigx2>1)=NaN; dBZ2(sigx2>1)=NaN;
yz2=dBZ2-132+26*sigx2.^2;
Ql2=(A*DVbar2).^(-3).*10.^(yz2/10);
Ql2(Ql2>10)=NaN;%g/m^3
Fq2=Ql2.*((DVbar2+B/A).*exp(-3*sigx2.^2)-B/A)*1e-6*1000*3600;%mm/hr
r02=1e3*(A*DVbar2+B).*exp(-6.5*sigx2.^2);%drizzle radius in mm

if newfigures; figure(4); dock; end; clf
plot(Ql.*detectm,ht,'linewidth',1.4);
hold on
plot(Fq,ht,'--','linewidth',1.4)
plot(r0,ht,'ro','linewidth',1.4)
ylabel('Height (m)');
xlabel({'Liquid Water (g/m^3)', 'Rainrate (mm/h)', 'Drizzle radius (mm)'});
axis([0 .25 0 2]);
% 2nd low-power noise threshold method 
plot(Ql2.*detectm,ht,'linewidth',1);
plot(Fq2,ht,'--','linewidth',1)
plot(r02,ht,'ro','linewidth',1)

