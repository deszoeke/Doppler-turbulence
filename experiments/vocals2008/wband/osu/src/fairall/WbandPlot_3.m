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
disp('flocmoi');

motionread=1;%set to zero to skip motion correction

ddds = input('Input yearday  ','s');
ddd = str2num(ddds);
shrs = input('Input start hr  ','s');
shr = str2num(shrs);
ehrs = input('Input end hr  ','s');
ehr = str2num(ehrs);
%eval(['cd ', way_raw_data_wband]);%eval(['cd ', 'C:\data\cwf\data\vocals\']);
save dew1 ddds ddd shrs ehrs shr ehr
%eval(['cd ', 'F:\VOCALS_2008\RHB\radar\wband\']);

%way_raw_images_wband= 'F:\VOCALS_2008\RHB\radar\wband\Raw_images\';
%way_raw_images_wband= 'c:\data\cwf\data\vocals\';

%ncload('C:\Data\VOCALS\Data\wband\NetCDFRadarmfiles\WbandMom.nc')  % First rain shower

%20082821346MMCRMom.nc  20082821400MMCRMom.nc% 2nd rain shower

ht = (1:1:120)*.025; % range gates in km
fn = [way_raw_data_wband '\20082821400MMCRMom.nc'];
fn(41:43) = ddds;

cldtop=[];
for jam = shr:ehr
    save dew2 jam
    %close all
    clear hh y z dt jd
    read_parameters;
    load dew1
    load dew2
    
    
    %way_raw_images_wband= 'F:\VOCALS_2008\RHB\radar\wband\Raw_images\';
    way_raw_images_wband= 'C:\Data\cwf\DATA\VOCALS_2008\cband_images';
    
    %ncload('C:\Data\VOCALS\Data\wband\NetCDFRadarmfiles\WbandMom.nc')  % First rain shower
    
    %20082821346MMCRMom.nc  20082821400MMCRMom.nc% 2nd rain shower
    
    ht = (1:1:120)*.025; % range gates in km
    %fn(28:30) = ddds;
    
    if jam<10,
        hr=['0' num2str(jam)];
    else
        hr=num2str(jam);
    end;
    fn = [way_raw_data_wband '\2008' ddds hr '00MMCRMom.nc'];
    FID = fopen(fn,'r');
    if FID > 0  % Test to see if file exists
        ncload(fn)
        
        tt = time_offset + base_time;
        % calculate time form time_offset and base_time
        for i = 1:length(tt)
            %X = base_time;
            y(i) = tt(i)/(60*60*24);
            z(i) = y(i)+datenum([1970,1,1,0,0,0]);
            dt = datestr(z(i), 31);    %2008-10-07 16:35:53
            yr = str2num(dt(1:4));
            month = str2num(dt(6:7));
            day = str2num(dt(9:10));
            hour = str2num(dt(12:13));
            minutes = str2num(dt(15:16));
            secs = str2num(dt(18:19));
            jd(i) = datenum(yr,month,day,hour,minutes,secs)-datenum(yr-1,12,31);
            %correct for
        end
        hh = (jd-(fix(jd(1))))*24; % subtract off JD and convert to an hour
        
        Ref = Reflectivity'; % Ref(height, time)
        SW = SpectralWidth';
        DVa = MeanDopplerVelocity';
        SNR = SignalToNoiseRatio';
        DV=[];
        if motionread%correct doppler for
            read_Wband_motion
            DV_mot=interp1(hhk,vkons(9,:),hh);
            ii=find(isnan(DV_mot));DV_mot(ii)=0;
            for im=1:length(ht);
                DV(im,:)=DVa(im,:)+DV_mot;
            end;
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
        ['Hour= ' num2str((hh(1))) '  Mean cloud top ht= ' num2str(httb) '  std cloud top height = ' num2str(sigh)]
        ['Mean dBZ= ' num2str(reftb) '   Mean dBZ-noisedBZ= ' num2str(reft2b) '   Mean Vel= ' num2str(wvelb)]
        cldtop=[cldtop' [hh(1) httb sigh (y(4)-y(3))/2 reftb reft2b wvelb nntb]']';
        
        figure(4)
        subplot(3,1,1),imagesc(hh,ht, Ref(1:120,:),[-60, 20]);  %
        axis xy
        ylabel('Ht (km)')
        tit = {[cruise year ' W-Band'' YD ' ddds ' Hr ' hr]};
        title(tit)
        h = colorbar;
        axes(h)
        ylabel('Ref dBZ')
        
        subplot(3,1,2),imagesc(hh,ht, DV(1:120,:),[-7, 7]);  %
        axis xy
        ylabel('Ht (km)')
        h = colorbar;
        axes(h)
        ylabel('DopVel (m/s)')
        
        
        subplot(3,1,3),imagesc(hh,ht, SW(1:120,:),[0, 5]);  %
        axis xy
        ylabel('Ht (km)')
        xlabel('Hrs (UTC)')
        h = colorbar;
        axes(h)
        ylabel('SpecWid (m/s)')
        print('-djpeg90 ',[way_raw_images_wband '\Wband_ref_' num2str(ddd), '_' hr '.jpg']);
        
        figure(5)
        clf
        subplot(3,1,1),imagesc(hh,ht, SNR(1:120,:),[-25, 20]);  %
        axis xy
        ylabel('Ht (km)')
        %axis([15 15.25 0 3])
        tit = {[cruise year ' W-Band'' YD ' ddds ' Hr ' hr]};
        title(tit)
        h = colorbar;
        axes(h)
        ylabel('SNR dBZ')
        
        subplot(3,1,2),imagesc(hh,ht, DV(1:120,:),[-3, 5]);  %
        axis xy
        ylabel('Ht (km)')
        %axis([15 15.25 0 3])
        h = colorbar;
        axes(h)
        ylabel('DopVel (m/s)')
        %
        %         subplot(3,1,3),imagesc(hh,ht, DVa(1:120,:),[-3, 5]);  %
        %         axis xy
        %         ylabel('Ht (km)')
        %         %axis([15 15.25 0 3])
        %         h = colorbar;
        %         axes(h)
        %         ylabel('DopVel (m/s)')
        
        subplot(3,1,3),imagesc(hh,ht, SW(1:120,:),[0, 5]);  %
        axis xy
        ylabel('Ht (km)')
        xlabel('Hrs (UTC)')
        h = colorbar;
        axes(h)
        ylabel('SpecWid (m/s)')
        print('-djpeg90 ',[way_raw_images_wband '\Wband_snr_' num2str(ddd), '_' hr '.jpg']);
        
        fclose(FID);
        
    end % test for file
    
end % hr loop

%%%%%%%%%%%%%%  Some diagnostics CWF added on the vocals cruise
ceilo=1;
if ceilo
    zz=load([way_proc_data_ceilo, cruise, year, 'ceilo_h_a.txt']);
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
end;

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
    
    save (['f:\VOCALS08_wband\reports\Radar_Cldtop_' ddds '.txt'],'cldtop','-tabs','-ascii')
end;
detect=zeros(size(Ref));
for i=1:length(Ref) %loop time i
    kk=find(Ref(:,i)+39-20*log10(ht')>8 & ht'>.1 & ht'<3 );
    detect(kk,i)=1;
end;
xRef=Ref.*detect;
xDV=DV.*detect;
xSW=SW.*detect;
xZZ=10.^(Ref/10).*detect;

dk=length(Ref)-1;kk=1;
xRefm=mean(xRef(:,kk:kk+dk)');
xZZm=mean(xZZ(:,kk:kk+dk)');
xDVm=mean(xDV(:,kk:kk+dk)');
xSWm=mean(xSW(:,kk:kk+dk)');
detectm=mean(detect(:,kk:kk+dk)');
mm=find(detectm<5/dk);xRevm(mm)=NaN;xDVm(mm)=NaN;xZZm(mm)=NaN; % xRevm-->xRefm?
dZBbar=10*log10(xZZm)./detectm;mm=find(dZBbar<-60);dZBbar(mm)=NaN;

ij=find(ht>2); % above cloud
wmm=mean(mean(DV(ij,kk:kk+dk)'));wdm=mean(mean(SW(ij,kk:kk+dk)'));
figure;
subplot(1,2,1);plot(xRefm./detectm,ht,mean(Ref(:,kk:kk+dk)'),ht,dZBbar,ht,'--',max(-60,-39+20*log10(ht')),ht,[-17 -17],[0 3],'--',[0 0],[0 3]);Ylabel('Height (m)');xlabel('dBZ');%axis([-60 20 0 2]);
subplot(1,2,2);plot(detectm,ht);xlabel('Return Fraction');%axis([0 1 0 2]);

% noise-filtered
DVbar=xDVm./detectm-wmm;
dBZ=xRefm./detectm;
Widthbar=xSWm./detectm;
wdbm=nanmean(Widthbar(ij));

% noisier, average weak Z returns too
DVbar2=mean(DV(:,kk:kk+dk)')-wmm;
dBZ2=mean(Ref(:,kk:kk+dk)');
ZZ=mean(10.^((Ref(:,kk:kk+dk)')/10));
dBZbar2=10*log10(ZZ);
Widthbar2=mean(SW(:,kk:kk+dk)');

figure;
subplot(1,2,1);
plot(xDVm./detectm-wmm,ht,mean(DV(:,kk:kk+dk)')-wmm,ht,'--',[0 0],[0 3]);Ylabel('Height (m)');xlabel('Mean Velocity (m/s)');%axis([-1 5 0 2]);
subplot(1,2,2);plot(xSWm./detectm,ht,mean(SW(:,kk:kk+dk)'),ht,'--',[0 0],[0 3]);xlabel('Doppler Width (m/s)');%axis([0 1 0 2]);

%From Frisch et al 1995
A=1.2e-4;
B=1e-5;

%sigx=sqrt(Widthbar.^2-wdbm^2)./DVbar;
sigx=sqrt(log(1+(Widthbar.^2-wdbm^2)./(DVbar+B/A).^2));
ii=find(sigx>1);sigx(ii)=NaN;DVbar(ii)=NaN;
yz=dBZ-132+26*sigx.^2;
Ql=(A*DVbar).^(-3).*10.^(yz/10);ii=find(Ql>10);Ql(ii)=NaN;%g/m^3
Fq=Ql.*((DVbar+B/A).*exp(-3*sigx.^2)-B/A)*1e-6*1000*3600;%mm/hr
r0=1e3*(A*DVbar+B).*exp(-6.5*sigx.^2);%drizzle radius in mm

% sigx2=sqrt(Widthbar2.^2-wdm^2)./DVbar2;
sigx2=sqrt(log(1+(Widthbar2.^2-wdm^2)./(DVbar2+B/A).^2));
ii=find(sigx2>1);sigx2(ii)=NaN;DVbar2(ii)=NaN;dBZ2(ii)=NaN;
yz2=dBZ2-132+26*sigx2.^2;
Ql2=(A*DVbar2).^(-3).*10.^(yz2/10);ii=find(Ql2>10);Ql2(ii)=NaN;%g/m^3
Fq2=Ql2.*((DVbar2+B/A).*exp(-3*sigx2.^2)-B/A)*1e-6*1000*3600;%mm/hr
r02=1e3*(A*DVbar2+B).*exp(-6.5*sigx2.^2);%drizzle radius in mm

figure;
subplot(1,2,1);
plot(Ql.*detectm,ht);Ylabel('Height (m)');xlabel('Liquid Water (g/m^3)');axis([0 .25 0 3]);
subplot(1,2,2);plot(Fq,ht);xlabel('Rainrate (mm/h)');axis([0 .3 0 3]);

