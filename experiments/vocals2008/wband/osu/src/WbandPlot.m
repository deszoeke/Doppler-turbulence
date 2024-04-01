% DEW KM read in W-Band nc moment data and plot
% range gates are 25 m
% timing is every 0.25 sec
close all
clear
%run('E:\VOCALS_2008\RHB\Scientific_analysis\programs\VOCALS2008_programs_cwf\read_parameters')
read_parameters

ls(way_raw_data_wband)
ddd=input('Input yearday  ');
ddds=sprintf('%03i',ddd);
shr=input('Input start hr  ');
ehr=input('Input end hr  ');
%eval(['cd ', 'H:\RHB\radar\wband\']);
%save dew1 ddds ddd shr ehr

%20082821346MMCRMom.nc  20082821400MMCRMom.nc% 2nd rain shower

%ht=(1:1:120)*.025; % range gates in km
cldtop=[];
suffix='MMCRMom.nc';
for hrn=shr:ehr
    hourstring=sprintf('%02i',hrn);
    fname=ls([way_raw_data_wband yyyy ddds hourstring '*' suffix]);
    if isempty(fname)
        disp(['No files found for hour ' hourstring])
        continue
    end
    count=0; % length of variables so far...
    for fi=1:size(fname,1) % loop reading all files in the hour (broken)
        fn=[way_raw_data_wband fname(fi,:)];
        if ~exist(fn,'file')
            disp([fname(fi,:) ' not found.']) % should never happen!
            continue
        end
        % file exists
        %ncload(fn) % broken on Simon's MATLAB R2007b
        % so switch to SNCTOOLS:
        % read scalars
        base_time=nc_varget(fn,'base_time'); % bug for multiple files in 1 hour!
        % read 1D vectors
        z=nc_varget(fn,'Heights');
        temp=nc_varget(fn,'time_offset')';
        num=length(temp);
        time_offset(count+1:count+num)=temp;
        tt(count+1:count+num)=base_time+temp;
        % read 2D variables
        Ref(count+1:count+num,:)=nc_varget(fn,'Reflectivity')+dB_offset; % nc_varget_lapxm?
        %+dB_offset 2013-04-28 SPdeS (1.72 dB more liberal cloud detection was used for the 2011 paper);
        SNR(count+1:count+num,:)=nc_varget(fn,'SignalToNoiseRatio');
        SW(count+1:count+num,:)=nc_varget(fn,'SpectralWidth');
        DV(count+1:count+num,:)=nc_varget(fn,'MeanDopplerVelocity');
        ncclose('all');
        count=count+num; % increment length of variables
    end % loop file(s) within hour
    
    % some transposes for now, but could be optimized not transposed

    ht=z(1,:)/1e3; % height in km (1st row)
    % calculate time from time_offset and base_time
    [yr,month,day,hour,minutes,secs]=datevec(tt/(60*60*24)+datenum([1970,1,1,0,0,0]));
    jd=datenum(yr,month,day,hour,minutes,secs)-datenum(yr,0,0); % decimal yearday of current year
    hh=(jd-(fix(jd(1))))*24; % subtract off JD and convert to decimal hour
    
    Ht=repmat(ht,[size(Ref,1) 1]); % height matrix
    kk=Ref+39-20*log10(Ht)>8 &...
        Ht>.5 &...
        Ht<3 &...
        Ref>-30; % logical
    nnt=sum(kk,2); % number of ranges satisfying kk condition per time
    whenk=nnt>0;   % logical times when there are kk selected
    [junk,klast]=max(kk(whenk,:),[],2);
    htt(whenk)=ht(klast); % the highest ht that satisfies kk conditon
    reft(whenk)=mean(10.^(Ref(kk(whenk,:))/10));
    reft2(whenk)=mean(10.^((Ref(kk(whenk,:))+39-20*log10(Ht(kk(whenk,:))))/10));
    wvel(whenk)=mean(DV(kk(whenk,:)));
    htt(~whenk)=NaN;
    reft(~whenk)=NaN;
    reft2(~whenk)=NaN;
    wvel(~whenk)=NaN;
%   ^ vectorized above ^
%     for k=1:length(Ref); % time index
%         kk=find(Ref(:,k)+39-20*log10(ht')>8 &...
%                 ht'>.5 &...
%                 ht'<3 &...
%                 Ref(:,k)>-30                    );
%         if ~isempty(kk)
%             htt(k)=ht(kk(length(kk)));
%             reft(k)=mean(10.^(Ref(kk,k)/10));
%             reft2(k)=mean(10.^((Ref(kk,k)+39-20*log10(ht(kk)'))/10));
%             wvel(k)=mean(DV(kk,k));
%             nnt(k)=length(kk);
%         else
%             htt(k)=NaN;
%             reft(k)=NaN;
%             reft2(k)=NaN;
%             wvel(k)=NaN;
%             nnt(k)=0;
%         end;
%     end;
    y=prctile(htt,[10 15 50 85 90]);
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

    figure(4)
    subplot(3,1,1),imagesc(hh,ht, Ref',[-60, 20]);  %
    axis xy
    ylabel('Ht (km)')
    title({[cruise year ' W-Band'' YD ' ddds ' Hr ' hourstring]})
    h=colorbar;
    set(get(h,'ylabel'),'string','Ref dBZ')

    subplot(3,1,2),imagesc(hh,ht, DV',[-7, 7]);  %
    axis xy
    ylabel('Ht (km)')
    h=colorbar;
    set(get(h,'ylabel'),'string','DopVel (m/s)')

    subplot(3,1,3),imagesc(hh,ht, SW',[0, 5]);  %
    axis xy
    ylabel('Ht (km)')
    xlabel('Hrs (UTC)')
    h=colorbar;
    set(get(h,'ylabel'),'string','SpecWid (m/s)')
    print('-djpeg90 ',[way_raw_images_wband '\Wband_ref_' num2str(ddd), '_' hourstring '.jpg']);

    figure(5)
    clf
    subplot(3,1,1),imagesc(hh,ht, SNR',[-25, 20]);  %
    axis xy
    ylabel('Ht (km)')
    %axis([15 15.25 0 3])
    title({[cruise year ' W-Band'' YD ' ddds ' Hr ' hourstring]})
    h=colorbar;
    set(get(h,'ylabel'),'string','SNR dBZ')

    subplot(3,1,2),imagesc(hh,ht, DV',[-3, 5]);  %
    axis xy
    ylabel('Ht (km)')
    %axis([15 15.25 0 3])
    h=colorbar;
    set(get(h,'ylabel'),'string','DopVel (m/s)')

    subplot(3,1,3),imagesc(hh,ht, SW',[0, 5]);  %
    axis xy
    ylabel('Ht (km)')
    xlabel('Hrs (UTC)')
    h=colorbar;
    set(get(h,'ylabel'),'string','SpecWid (m/s)')
    print('-djpeg90 ',[way_raw_images_wband '\Wband_snr_' num2str(ddd), '_' hourstring '.jpg']);

end % loop over hours

%%%%%%%%%%%%%%  Some diagnostics CWF added on the vocals cruise
ceilo=1;
if length(cldtop)>22;
    %ddd=327;ddds=num2str(ddd);cldtop=load(['cldtop_' ddds '.txt']);
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

    jdt=ddd+cldtop(:,1)/24+.5/24;cd
    cldtop(cldtop(:,8)<500,2:7)=NaN;
    figure;
    subplot(2,1,1);plot(jdt,cldtop(:,2),'-o',jdt,cldtop(:,2)+cldtop(:,3),'--',jdt,cldtop(:,2)-cldtop(:,3),'--',jdt,cldtop(:,2)-cldtop(:,8)/12600*25/1e3,'x-');
    if ceilo;hold on;plot(jdh,cloudhgt/1000,'-d',jdh,hgts1/1000,'.',jdh,hgts2/1000,'.');axis([ddd ddd+1 .5 1.8]);end;
    Ylabel('Cloud top/base height (km)');
    subplot(2,1,2);plot(jdt,cldtop(:,5),'-o',jdt,10*cldtop(:,7),'-x',[ddd ddd+1],[0 0],'-' ,[ddd ddd+1],[-17 -17],'--' );
    Ylabel('dBZ and 10*W (m/s)');
    xlabel('Julian Day (2008)');
    
    figure;
    plot(jdt,cldtop(:,2)-cloudhgt(floor(jdh)==ddd)/1000);axis([ddd ddd+1 0 .5]);
    Ylabel('Cloud Thickness (m)');
    xlabel('Julian Day (2008)');

    save (['c:\data\cwf\doc\climate\vocals\field_program\RHB\reports\' ddds '\cldtop_' ddds '.txt'],'cldtop','-tabs','-ascii')
end;

detect=Ref+39-20*log10(Ht)>8 & Ht>.1 & Ht<3; % logical
kk=find(detect);
% ^ vectorized ^
%  detect=zeros(size(Ref));
%  for i=1:length(Ref)
%       kk=find(Ref+39-20*log10(Ht')>8 & Ht'>.1 & Ht'<3 );
%       detect(kk,i)=1;
%  end;
 xRef=Ref.*detect;
 xDV=DV.*detect;
 xSW=SW.*detect;
 xZZ=10.^(Ref/10).*detect;

 % (hourly) time means for each range: (time dim. already 1)
 dk=length(Ref)-1;
 % kk=1; % not used
 xRefm=mean(xRef,1);
 xZZm=mean(xZZ,1);
 xDVm=mean(xDV,1);
 xSWm=mean(xSW,1);
 detectm=mean(detect,1); % fraction of times something was detected
 
 mm=find(detectm<5/dk);
 xRevm(mm)=NaN;
 xDVm(mm)=NaN;
 xZZm(mm)=NaN;
 dZBbar=10*log10(xZZm)./detectm;
 dZBbar(dZBbar<-60)=NaN;
  
 figure;
 subplot(1,2,1);
 plot(xRefm./detectm,ht,mean(Ref),ht,dZBbar,ht,'--',max(-60,-39+20*log10(ht')),ht,[-17 -17],[0 3],'--',[0 0],[0 3]);
 ylabel('Height (m)');
 xlabel('dBZ');%axis([-60 20 0 2]);
 subplot(1,2,2);
 plot(detectm,ht);
 xlabel('Return Fraction');%axis([0 1 0 2]);
 
 wmm=mean(mean(DV(ht>2,:)));
 wdm=mean(mean(SW(ht>2,:)));
 DVbar=xDVm./detectm-wmm;
 dBZ=xRefm./detectm;
 Widthbar=xSWm./detectm;
 wdbm=nanmean(Widthbar(ht>2));
 
 DVbar2=mean(DV)-wmm;
 dBZ2=mean(Ref);
 ZZ=mean(10.^(Ref/10));
 dBZbar2=10*log10(ZZ);
 Widthbar2=mean(SW);
 
 figure;
 subplot(1,2,1);
 plot(xDVm./detectm-wmm,ht,mean(DV)-wmm,ht,'--',[0 0],[0 3]);
 ylabel('Height (m)');
 xlabel('Mean Velocity (m/s)');%axis([-1 5 0 2]);
 subplot(1,2,2);
 plot(xSWm./detectm,ht,mean(SW),ht,'--',[0 0],[0 3]);
 xlabel('Doppler Width (m/s)');%axis([0 1 0 2]);
 
 %From Frisch et al.- 1995
 A=1.2e-4;
 B=1e-5;
 
 %sigx=sqrt(Widthbar.^2-wdbm^2)./DVbar;
 sigx=sqrt(log(1+(Widthbar.^2-wdbm^2)./(DVbar+B/A).^2));
 %                             ^^^ not in eqn. 14 Frisch et al. 1995
 sigx(sigx>1)=NaN;
 DVbar(sigx>1)=NaN;
 yz=dBZ-132+26*sigx.^2;
 Ql=(A*DVbar).^(-3).*10.^(yz/10);
 Ql(Ql>10)=NaN;%g/m^3
 Fq=Ql.*((DVbar+B/A).*exp(-3*sigx.^2)-B/A)*1e-6*1000*3600;%mm/hr
 r0=1e3*(A*DVbar+B).*exp(-6.5*sigx.^2);%drizzle radius in mm
 
% sigx2=sqrt(Widthbar2.^2-wdm^2)./DVbar2;
 sigx2=sqrt(log(1+(Widthbar2.^2-wdm^2)./(DVbar2+B/A).^2));
 sigx2(sigx2>1)=NaN;
 DVbar2(sigx2>1)=NaN;
 dBZ2(sigx2>1)=NaN;
 yz2=dBZ2-132+26*sigx2.^2;
 Ql2=(A*DVbar2).^(-3).*10.^(yz2/10);
 Ql2(Ql2>10)=NaN;%g/m^3
 Fq2=Ql2.*((DVbar2+B/A).*exp(-3*sigx2.^2)-B/A)*1e-6*1000*3600;%mm/hr
 r02=1e3*(A*DVbar2+B).*exp(-6.5*sigx2.^2);%drizzle radius in mm
 
 figure;
 subplot(1,2,1);
 plot(Ql.*detectm,ht);
 ylabel('Height (m)');
 xlabel('Liquid Water (g/m^3)');
 axis([0 .25 0 3]);
 subplot(1,2,2);
 plot(Fq,ht);
 xlabel('Rainrate (mm/h)');
 axis([0 .3 0 3]);
