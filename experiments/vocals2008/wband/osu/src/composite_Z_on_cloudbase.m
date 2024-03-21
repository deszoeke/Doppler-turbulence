% Composite W-band reflectivity Z on column max. reflectivity CMR
read_parameters;
load([way_proc_data_wband '1min_stat/Z_1min.mat']);

% Define dimensions
% vertical coordinate (center of ranges, m)
% h=load([way_proc_data_wband '1min_stat/height.txt']);
% % dBZ and w bin edges
% Zbins=[-inf -40:2:30 inf];
% wbins=[-inf -5:0.2:5 inf];

nheight=length(Z.height); % 75
nmin=60;
nhour=24;
nZbin=length(Z.bins);

% load 1-min W-band cloud top height (Simon de Szoeke, Oregon State U.)
% base_time=1225916508;
A=load([way_proc_data_wband 'cloudheight/CloudHeight_1min_2008310-336.txt']);
A(A==-999)=NaN; % NaN out missing value of -999
% already synchronzied to time_yday from 1min_stat/Z_1min.mat
% time=A(:,1); % seconds since base_time
cloud_fraction=A(:,2);
cloud_top_mean=A(:,3);
cloud_top_std=A(:,4);
cloud_top_val15=A(:,5);
cloud_top_median=A(:,6);
cloud_top_val85=A(:,7);
cloud_top_valmean=A(:,8); %mean of tops within 15-85 percentile
num_returns=A(:,9);
num_cloud_returns=A(:,10);
num_clear_returns=A(:,11);
num_nonmet_returns=A(:,12);
flag=A(:,13);
% flag values
% 0 no errors detected, cloud top computed.
% 1 cloud fraction=0, cloud top set to NaN.
% 2 cloud fraction>=0 but found no 3 cloud returns with consistent height,
%   possible glitch, provisional cloud top computed.
% NaN out percentile-based cloud tops=0 using flag
lf=logical(flag);
cloud_top_val15(lf)=NaN;
cloud_top_val85(lf)=NaN;
cloud_top_valmean(lf)=NaN;
cloud_top_median(lf)=NaN;
% time_yday=datenum(0,0,0,0,0,base_time+time_aggreg)+datenum(1970,1,1,0,0,0)-datenum(2008,1,0,0,0,0); %yearday

% load 10-min cloud boundaries
ncpath='/Users/sdeszoek/Data/cruises/VOCALS_2008/RHB/radar/wband/Processed/cloudheight/'; % note netcdflib does not accept ~/ as a path
ncf=[ncpath 'VOCALS2008CloudBoundaries10min_0_2.nc'];
yday=nc_varget(ncf,'yday');
lat=nc_varget(ncf,'lat');
lon=nc_varget(ncf,'lon');
cloudtop=nc_varget(ncf,'cloudtop');
sondecloudtop=nc_varget(ncf,'sondecloudtop');
sondeinterpflag=nc_varget(ncf,'sondeinterpflag');
cloudbase=nc_varget(ncf,'cloudbase');
ncloud=nc_varget(ncf,'numceilcloud');
time10yday=yday;

% load 10-min CL31 cloud base height and cloud fraction
A=load('/Users/sdeszoek/Data/cruises/VOCALS_2008/RHB/ceilometer/Processed/VOCALS2008ceilo_10min_a.txt');
clb10.yday=A(:,1);
clb10.cloudfrac=A(:,9);
clb10.height=A(:,13); % 85 percentile 1st cloud base height
clb10.height(clb10.height>2000)=NaN; % filter above 2000m

% load 10-min LWP (Paquita Zuidema, U. Miami)
lwpfile=[way_proc_data_umiami 'vocals_rhb_lwp10min_v1.nc'];
lwp.yday =nc_varget(lwpfile,'jd_jdfrac');
lwp.lwp  =nc_varget(lwpfile,'lwp'); % g/m^2, 2-band mailbox physical retrieval (Zuidema 2005)
lwp.nobs =nc_varget(lwpfile,'nlwp');
lwp.adlwp=nc_varget(lwpfile,'adlwp');
% QC data
lwp.lwp(lwp.nobs<10)=NaN;
lwp.adlwp(lwp.adlwp>1e10 | lwp.adlwp<0)=NaN;

% load flux file for position
A=load([way_proc_data_flux 'cat/flux_5hf_VOCALS2008.txt']);
flux.yday=A(:,1);
flux.lat=A(:,17);
flux.lon=A(:,18);

% load sondes
load([way_raw_data_balloon(1:end-4) 'mat/wtec_edt_interp2z.mat'])
yday_sonde=datenum(year,month,day,hour,minute,0)-datenum('0-jan-2008');

% SYNCHRONIZE 10-min data to time10yday...
% synchronize LWP and Ztopmax10(time10yday)
[ii,loc]=ismember(round(lwp.yday*60*24/10),round(time10yday*60*24/10));
% time10yday(loc(ii)) approx equals lwp.yday(ii)
lwp10=NaN+time10yday;
adlwp10=lwp10;
lwp10(loc(ii))=lwp.lwp(ii);
adlwp10(loc(ii))=lwp.adlwp(ii);

% define a noise level for the reflectivity based on Ken Moran's manual
radar_const=19.613; % dB_
dB_offset=-1.72; % dB
min_detectable_signal=-113.304+dB_offset; % dBmW
analog_noise=-120.74;  % peak, dBmW
digital_noise=-115.4;   % peak, dBmW
threshold_margin=3.5; % dB
noisefloor=20*log10(Z.height)+radar_const+digital_noise+threshold_margin;
% do not plot colors/contours below the noise threshold
isnoise=Z.mean<repmat(noisefloor',[size(Z.mean,1) 1]);
Z.plotmean=Z.mean;
Z.plotmean(isnoise)=NaN;
Z.plotvar=Z.var;
Z.plotvar(isnoise)=NaN;

% Column maximum reflectivity (CMR)
% vertical max reflectivity of 1-minute means
[zmx,k]=max(Z.mean,[],2);

% composite 1-min Z profiles on 10-minute cloudbase
basebin=[600 800 1000 1200 1400 1600 1800]; % 6 bins defined between these 7 endpoints
ntime=length(time_yday);
nheight=length(Z.height);
for indbin=1:length(basebin)-1
    ii=cloudbase>=basebin(indbin) & cloudbase<basebin(indbin+1);
    jj=false(size(Z.time_yday));
    fit=find(ii);
    for ti=1:sum(ii)
        jj=jj + Z.time_yday>=yday(fit(ti)) & Z.time_yday<yday(fit(ti)+1);
    end
    countinbin(indbin)=sum(jj);
    Z.basebin(indbin,:)=nanmean(Z.mean(jj,:)); % simple level-mean
    Z.basecfad(indbin,:,:)=nansum(Z.cfad(jj,:,:));
end

h=area([-60; noisefloor],Z.height([1 1:end]),2000,'facecolor',0.7+[0 0 0],'edgecolor','none');
hold on
plot(Z.basebin,Z.height,'linewidth',2)
axis([-50 20 0 2e3])
set(gca,'fontsize',16)
ylabel('altitude (m)')
xlabel('cloud base height composite mean reflectivity (dBZ)')
%print('-dpng',[way_proc_images_wband 'composite_Z_on_cloudbase.png'])

% CFAD of Z
for indbin=1:6
    ax(indbin)=subplot(2,3,indbin,'align');
    
    [cc,hc]=contourf(Z.bins(2:end-1)+(Z.bins(3)-Z.bins(2))/2,Z.height/1e3,squeeze(Z.basecfad(indbin,:,2:end-1)),40);
    set(hc,'edgecolor','none')
    colormap(b2rcolormap(40));
    hb=colorbar('southoutside');
    set(get(hb,'xlabel'),'string','count')
    title(sprintf('z_{base}=[%3.1f %3.1f) km',basebin(indbin+[0 1])/1e3));
    hold on
    plot(noisefloor,Z.height/1e3,'w-')
    plot(Z.basebin(indbin,:),Z.height/1e3,'m-')
    xlabel('dBZ')
    ylabel('height (km)')
end
%print('-dpng',[way_proc_images_wband 'comp_Zcfad_on_cloudbase.png'])

% movie of CFADs
edges=400:100:1700;
count=histc(cloudbase,edges);

hf=figure; dock(-hf);
clf
ax(1)=subplot(2,1,1); hold on;
    hb=colorbar('southoutside','fontsize',14);
    set(hb,'position',get(hb,'position')-[0 .04 0 0])
    set(get(hb,'xlabel'),'string','count','fontsize',16)
    hold on; delete(hb);
ax(2)=subplot(3,1,3);
hbar=bar(edges/1e3,count/60,'histc');
set(hbar,'facecolor',0.7*[1 1 1])
axis([.3 1.8 0 20])
set(gca,'tickdir','out','fontsize',14)
ylabel('frequency (hours)')
xlabel('cloud base height (km)') 
set(gca,'tickdir','out')
hold on

indbin=1;
    axes(ax(2));
    ha=plot(edges(indbin)+25,-3,'k^','markersize',8,'markerfacecolor','k','clipping','off');
    
    axes(ax(1));
    ii=cloudbase>=edges(indbin) & cloudbase<edges(indbin+1);
    cfad(indbin,:,:)=nansum(Z.cfad(ii,:,2:end-1),1);
    if any(any(cfad(indbin,:,:)~=0))
        [cc,hc]=contourf(Z.bins(2:end-1)+(Z.bins(3)-Z.bins(2))/2,Z.height/1e3,squeeze(cfad(indbin,:,:)),40);
        set(hc,'edgecolor','none')
        ht=text(1,1,'');
    else
        hc=image(Z.bins(2:end-1)+(Z.bins(3)-Z.bins(2))/2,Z.height/1e3,squeeze(cfad(indbin,:,:))); % blue screen
        axis xy; axis([-39.0000   31.0000    0.1336    1.9818]);
        ht=text(-4,1,'no scenes','fontsize',16,'color','w','horizontalalignment','center');
    end
    set(gca,'fontsize',16)
    colormap(b2rcolormap(40));
    title(sprintf('z_{base}=[%3d %d) m',edges(indbin+[0 1])));
    xlabel('dBZ')
    ylabel('height (km)')
    hb=colorbar('southoutside','fontsize',14);
    set(hb,'position',get(hb,'position')-[0 .04 0 0])
    set(get(hb,'xlabel'),'string','count','fontsize',16)
    hold on
    plot(noisefloor,Z.height/1e3,'w-')
    % make adjustments
    print('-djpeg',[way_proc_images_wband sprintf('CFAD_movie/ZCFAD%02dcloudbase%+03d.jpeg',indbin,edges(indbin))])
    frame(indbin)=getframe(hf);
    delete(hc,ha,ht,hb);

for indbin=2:length(edges)-1
    axes(ax(2));
    ha=plot(edges(indbin)+25,-3,'k^','markersize',8,'markerfacecolor','k','clipping','off');
    
    axes(ax(1));
    ii=cloudbase>=edges(indbin) & cloudbase<edges(indbin+1);
    cfad(indbin,:,:)=nansum(Z.cfad(ii,:,2:end-1),1);
    if any(any(cfad(indbin,:,:)~=0))
        [cc,hc]=contourf(Z.bins(2:end-1)+(Z.bins(3)-Z.bins(2))/2,Z.height/1e3,squeeze(cfad(indbin,:,:)),40);
        set(hc,'edgecolor','none')
        ht=text(1,1,'');
    else
        hc=image(Z.bins(2:end-1)+(Z.bins(3)-Z.bins(2))/2,Z.height/1e3,squeeze(cfad(indbin,:,:))); % blue screen
        axis xy; axis([-39.0000   31.0000    0.1336    1.9818]);
        ht=text(-4,1,'no scenes','fontsize',16,'color','w','horizontalalignment','center');
    end
    set(gca,'fontsize',16)
    colormap(b2rcolormap(40));
    title(sprintf('z_{base}=[%d %d]',edges(indbin+[0 1])));
%     xlabel('dBZ')
%     ylabel('height (km)')
    hb=colorbar('southoutside','fontsize',14);
    set(hb,'position',get(hb,'position')-[0 .04 0 0])
    set(get(hb,'xlabel'),'string','count','fontsize',16)
    hold on
    plot(noisefloor,Z.height/1e3,'w-')
    print('-djpeg',[way_proc_images_wband sprintf('cloudbase_movie/ZCFAD%02dcloudbase%+03d.jpeg',indbin,edges(indbin))])
    frame(indbin)=getframe(hf);
    delete(hc,ha,ht,hb);
end

movie2avi(frame,[way_proc_images_wband 'cloudbase_movie/ZCFAD_on_cloudbase.avi'],'fps',5);