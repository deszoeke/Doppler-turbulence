% plot_wband_nov20.m
% 2010-04-08 :: VOCALS 2008 :: Simon de Szoeke
%
% Plot the wband moments Z and w for Nov 20, 11-13 UTC
%

%cd /Users/sdeszoek/Data/cruises/VOCALS_2008/RHB/Scientific_analysis/programs/VOCALS2008_programs/wband
%addpath('~/Data/cruises/VOCALS_2008/RHB/Scientific_analysis/programs/VOCALS2008_programs/');
read_parameters;

% load flux file for position
A=load([way_proc_data_flux 'cat/flux_5hf_VOCALS2008.txt']);
flux.yday=A(:,1);
flux.lat=A(:,17);
flux.lon=A(:,18);

% All radar data files
momentfile=[way_raw_data_wband '20083251000MMCRMom.nc';
            way_raw_data_wband '20083251100MMCRMom.nc';
            way_raw_data_wband '20083251200MMCRMom.nc';
            way_raw_data_wband '20083251300MMCRMom.nc' ];

h=nc_varget_lapxm(momentfile(1,:),'Heights',[0 0],[1 -1])';

% After best recalibration, Ken Moran says to subtract 1.72 dB to all recorded values.
% reviever gain is 1.72 higher than previously estimated, which reduces dBZ
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

% loop files
%starter=find(base_time_yday>=310,1,'first');
countt=0;
for fi=1:size(momentfile,1) % fi>=2;
    %%prflname=[way_raw_data_wband momentfile(fi-1).name]; % previous file
    filename=momentfile(fi,:);
    
    year=2008;
    % found no need to read the end of the previous file to get an integral minute
    % base_time in seconds since 1970-1-1 00:00:00
    base_time=nc_varget(filename,'base_time');
    % time_offset in seconds
    time_offset=nc_varget(filename,'time_offset');
    nt=length(time_offset);
    time(countt+(1:nt))=int32(time_offset)+base_time;
    % base_time_mld in matlab datenumber
    base_time_mld=double(base_time)/86400 + datenum(1970,1,1,0,0,0);
    base_time_yday=base_time_mld-datenum(year,1,0,0,0,0);
    time_mld=double(time)/86400 + datenum(1970,1,1,0,0,0);
    time_yday=time_mld-datenum(year,1,0,0,0,0);
        
    % reflectivity matrix
    Z(countt+(1:nt),:)=nc_varget_lapxm(filename,'Reflectivity',[0 0],[-1 length(h)]);
    vel(countt+(1:nt),:)=nc_varget_lapxm(filename,'MeanDopplerVelocity',[0 0],[-1 length(h)]);
    width(countt+(1:nt),:)=nc_varget_lapxm(filename,'SpectralWidth',[0 0],[-1 length(h)]);
    
    % Read Kongsberg motion compensation file
    DV=datevec(base_time_mld);
    yyyy=sprintf('%04d',year);
    ddd= sprintf('%03i',floor(base_time_yday));
%     hh=  sprintf('%02d',floor(mod(base_time_yday,1)*24)); % subject to rounding error!
    DV=datevec(base_time_mld);
    hh=  sprintf('%02d',DV(4));
    kongfile=dir([way_raw_data_wband 'motion_adjT/' yyyy ddd hh '*Kongsberg_adjT.txt']);
    if isempty(kongfile)
        kongflag=1;
        fkongerr=fopen([way_raw_data_wband 'motion_adjT/kongsberg_log.txt'],'a+');
        fprintf(fkongerr,['No Kongsberg file in hour:' ddd ' ' hh]);
        fclose(fkongerr);
        H=repmat(h,[nt 1]); % no angle correction
    elseif length(kongfile)>1
        kongflag=2;
        fkongerr=fopen([way_raw_data_wband 'motion_adjT/kongsberg_log.txt'],'a+');
        fprintf(fkongerr,['Multiple files found: ' kongfile(:).name]);
        fclose(fkongerr);
        H=repmat(h,[nt 1]); % no angle correction
    else
        kongflag=0;
        kongfilename=[way_raw_data_wband 'motion_adjT/' kongfile(1).name];
        [kongtime,kongw,kongpitch,kongroll]=read_kongsberg(kongfilename);
        
        % synchronize Kongsberg and radar time (yearday)
        [junk,ii,jj]=unique(kongtime);
        kongw4radar(countt+(1:nt),1)=interp1(kongtime(ii),kongw(ii)    ,base_time_yday+time_offset/86400);
        pitch(countt+(1:nt),1)      =interp1(kongtime(ii),kongpitch(ii),base_time_yday+time_offset/86400); %degrees
        roll(countt+(1:nt),1)       =interp1(kongtime(ii),kongroll(ii) ,base_time_yday+time_offset/86400); %      
        % pitch/roll housekeeping vertical coordinate
        quad=sin(pitch/180*pi).^2+sin(roll/180*pi).^2; % equals sin(theta)^2
        %sintheta=sqrt(quad);
        %theta=asin(sintheta); % zenith angle (radians)
        costheta=sqrt(1-quad);
        if sum(isfinite(costheta))<max(3.5*60,length(costheta)/2)
            kongflag=3; % insufficient motion data
            H=repmat(h,[nt 1]); % no angle correction
            H(isfinite(costheta),:)=costheta(isfinite(costheta))*h; % correct where available
        else
            H=costheta*h'; % antenna zenith angle cosine correction
        end
    end
    countt=countt+nt;

end

% ship-motion-compensated vertical velocity +towards
w=vel+repmat(kongw4radar,[1,length(h)]);

mask=Z>=repmat(noisefloor',[length(Z) 1]);
Z(~mask)=NaN;
w(~mask)=NaN;
width(~mask)=NaN;
Z(:,h<170)=NaN;
w(:,h<170)=NaN;
width(:,h<170)=NaN;

figure(1); dock; clf
drizzle= [ 0.3137    0.3176    0.3137
           0.4092    0.4118    0.4092
           0.5046    0.5059    0.5046
           0.6000    0.6000    0.6000
           0.6637    0.6608    0.6824
           0.7275    0.7216    0.7647
           0.7912    0.7824    0.8471
           0.8549    0.8431    0.9294
           0.6745    0.6059    0.9647
           0.4941    0.3686    1.0000
           0.2706    0.3608    1.0000
           0.2000    0.6000    1.0000
           0         1.0000         0
           1.0000    1.0000         0
           1.0000    0.8000         0
           1.0000    0.6000         0 ];
colormap([drizzle; b2rcolormap(17)]);

ax(1)=subplot(2,1,1,'align');
hpc(1)=pcolor(rem(time_mld,1)*24,h/1e3,Z'); shading flat
axis xy
axis([10.5 14 0.1 2])
caxis([-50 110])
hb(1)=colorbar; set(hb(1),'ylim',[-50 30])
set(gca,'color','k','fontsize',14,'tickdir','out')
title('reflectivity (dBZ)')

ax(2)=subplot(2,1,2,'align');
hpc(2)=pcolor(rem(time_mld,1)*24,h/1e3,-w'); shading flat
axis xy
axis([10.5 14 0.1 2])
caxis([-18 6])
hb(2)=colorbar; set(hb(2),'ylim',[-6 6])
set(gca,'color','k','fontsize',14,'tickdir','out')
title('vertical velocity (m/s)')
xlabel('Nov 20 hour')
ylabel('height (km)')
set(gcf,'inverthardcopy','off','color','w')
% print('-dpng',[way_proc_images_wband 'wband_nov20.png'])
% delete(hpc(:))
% print('-depsc',[way_proc_images_wband 'wband_nov20.eps'])

% plot a period where the width becomes very large at cloud base
colormap([[0 0 0]; b2rcolormap(17)]);
ti=rem(time_mld,1)*24>=11.2 & rem(time_mld,1)*24<=11.6;
ax(1)=subplot(2,1,1,'align');
hpc(1)=imagesc((rem(time_mld(ti),1)*24-11)*60,h/1e3,-w(ti,:)'); shading flat
axis xy
axis([18 36 0.1 2])
caxis([-5.5 3])
hb(2)=colorbar; set(hb(2),'ylim',[-5.5 1],'fontsize',14)
set(gca,'color','k','fontsize',14,'tickdir','out')
title('vertical velocity (m/s)')
% xlabel('Nov 20 hour')
ylabel('height (km)')
set(gcf,'inverthardcopy','off','color','w')

ax(2)=subplot(2,1,2,'align');
hpc(2)=imagesc((rem(time_mld(ti),1)*24-11)*60,h/1e3,width(ti,:)');
axis xy
axis([18 36 0.1 2])
caxis([-0.25 4])
hb(2)=colorbar; set(hb(2),'ylim',[-0.25 4],'fontsize',14)
set(gca,'color','k','fontsize',14,'tickdir','out')
title('Doppler width (m/s)')
xlabel('Nov 20 11 UTC minute')
ylabel('height (km)')
set(gcf,'inverthardcopy','off','color','w')
% print('-depsc',[way_proc_images_wband 'vel_width_nov20.eps'])

% velocity vs. Z joint histograms
zedge=-40:2:20; % 30 bins
vedge=-3:0.2:6; % 45 bins;
count=zeros(120,45,30);
for i=1:30 % reflectivity
    ii=Z>=zedge(i) & Z<zedge(i+1);
    for j=1:45 % velocity
        jj=w>=vedge(j) & w<vedge(j+1);
        count(:,j,i)=sum(ii&jj);        
    end
end
pcolor(zedge,vedge,squeeze(sum(count(:,[1:end 1],[1:end 1]))))

%% parameterize velocity by dBZ after Frisch lognormal dist
% Z is dBZ
dBZ=-40:20;
sigx=0.35;
a=1.2e-4; % s
b=1.0e-5; % m
N=1.5e3; % drops per m^3

const=60*log10(2) + 10*log10(N) + 180*sigx^2/log(10) + 180;
r0=10.^((dBZ-const)/60);
Vd=(exp(13*sigx^2/2)*r0 - b)/a;
plot(dBZ,Vd,'r','linewidth',1)

%% plot histograms for different heights
hedge=[200:200:2000]; % 10 edges -> 9 bins
clf
for i=1:9 % height
    ax(i)=subplot(3,3,10-i);
    ii=h>=hedge(i) & h<hedge(i+1);
    contourf(zedge(1:end-1)+1,vedge(1:end-1)+0.1,squeeze(sum(count(ii,:,:))),20); shading flat;
    title(sprintf('[%03.1f %03.1f) km',hedge(i+[0 1])/1e3),'fontsize',14)
    str=sprintf('%2.1e',sum(sum(count(ii,:))));
    text(-10,-1.1,str)
    hold on
    bigsum=squeeze(sum(sum(count(ii,:,:)),3));
    plot(bigsum/max(bigsum)*50-40,vedge(1:end-1)+.1)
    plot(dBZ,Vd,'r','linewidth',1)
end
set(ax(:),'ylim',[-1.5 4.5])
set(get(ax(3),'xlabel'),'string','dBZ','fontsize',14)
set(get(ax(3),'ylabel'),'string','m s^{-1}','fontsize',14)
b=colormap(1-bone);
colormap([[1 1 1]; b(2:end,:)]);
orient tall
%print('-depsc',[way_proc_images_wband 'joint_w_Z_nov20.eps'])
%print('-dpng',[way_proc_images_wband 'joint_w_Z_nov20.png'])

for i=1:9 % height
    ax(i)=subplot(3,3,10-i);
    ii=h>=hedge(i) & h<hedge(i+1);
    plot(Z(:,ii),w(:,ii),'k.','markersize',1);
    title(sprintf('[%03.1f %03.1f) km',hedge(i+[0 1])/1e3))
end
set(get(ax(7),'xlabel'),'string','dBZ')
set(get(ax(7),'ylabel'),'string','m/s')
set(ax(:),'ylim',[-6 6],'xlim',[-40 20])

% CFAD
zedges=-50:1:30;
ntime=length(time_yday);
nheight=size(Z,2);
zcfad=zeros(length(zedges)-1,nheight);
for indbin=1:length(zedges)-1;
    ii=Z>=zedges(indbin) & Z<zedges(indbin+1);
    zcfad(indbin,:)=nansum(ii,1);
end
wedges=-6:.25:6;
wcfad=zeros(length(wedges)-1,nheight);
for indbin=1:length(wedges)-1;
    ii=w>=wedges(indbin) & w<wedges(indbin+1);
    wcfad(indbin,:)=nansum(ii,1);
end
wiedges=0:.075:6;
wicfad=zeros(length(wiedges)-1,nheight);
for indbin=1:length(wiedges)-1;
    ii=width>=wiedges(indbin) & width<wiedges(indbin+1);
    wicfad(indbin,:)=nansum(ii,1);
end

% also plot from plot_wband_nov18.m
b2rcolormap(24);
subplot(3,2,2,'align')
pcolor(zedges(1:end-1),h,log(zcfad'))
shading flat; hold on
plot(noisefloor,h,'y-')
set(gca,'ylim',[0 2000])
title({'2008 Nov 20' 'reflectivity log-CFAD'})
xlabel('dBZ')
subplot(3,2,4,'align')
pcolor(wedges(1:end-1),h,log(wcfad'))
shading flat; hold on
plot([0 0],[0 3000],'y')
set(gca,'ylim',[0 2000])
title('vertical velocity log-CFAD')
subplot(3,2,6,'align')
pcolor(wiedges(1:end-1),h,log(wicfad'))
shading flat; hold on
set(gca,'ylim',[0 2000],'xlim',[0 5])
title('Doppler vel. width log-CFAD')
xlabel('m s^{-1}')

set(gcf,'color','w','inverthardcopy','off')
orient tall
% print('-depsc',[way_proc_images_wband 'zwwidcfad_nov18_20.eps'])
% print('-dpng',[way_proc_images_wband 'zwwidcfad_nov18_20.png'])
% % delete expensive images
% print('-depsc',[way_proc_images_wband 'zwwidcfad_nov18_20_blank.eps'])
