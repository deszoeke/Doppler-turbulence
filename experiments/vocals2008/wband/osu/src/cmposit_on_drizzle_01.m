% composite Z and w from W-band radar on drizzle index.
% Iterate on drizzle index, shaped by case studies from Nicholas Elmer.
%
% VOCALS 2008 :: 2011-07-05 :: Simon de Szoeke

read_parameters;
put_gordon=[way_proc_images_wband 'Gordon2011/'];

load([way_proc_data_wband '1min_stat/Z_1min.mat']);
load([way_proc_data_wband '1min_stat/w_1min.mat']);
load([way_proc_data_wband '1min_stat/Dw_1min.mat']);

% Is CMR at cloud base?
ncpath='/Users/sdeszoek/Data/cruises/VOCALS_2008/RHB/radar/wband/Processed/cloudheight/';
% note netcdflib does not accept ~/ as a path
ver=1.0;
ncf=[ncpath 'VOCALS2008CloudBoundaries10min_v' sprintf('%3.1f',ver) '.nc'];
cldbnd_yday=nc_varget(ncf,'yday');
cloudtop=nc_varget(ncf,'cloudtop');
cloudbase=nc_varget(ncf,'cloudbase');

% synchronize to 10 minute
t10=floor(Z.time_yday*144)/144;
z10=NaN+t10;
for i=find(cldbnd_yday-315<1/144):length(cldbnd_yday)
    ii=abs(t10-cldbnd_yday(i))<1/288;
    if ~isempty(ii)
        z10(ii)=cloudtop(i);
    end
end

% height of col max refl = Z.height(k)
[zmx,k]=max(Z.mean,[],2);
zht=Z.height(k);
zht(zmx<-34)=NaN;

% plot(cldbnd_yday,cloudbase)
% hold on
% plot(cldbnd_yday,cloudtop)
% plot(Z.time_yday,Z.height(k),'r.','markersize',2)

plot(zht-z10,'.','markersize',1)

edges=-1800:25:100;
hst=histc(zht-z10,edges);

b=barh(edges/1e3,hst/max(hst),'histc');
hold on
set(b,'facecolor',[0.7 0.7 0.7])
plot(cumsum(hst)/sum(hst),[edges(2:end) 200]/1e3)
set(gca,'tickdir','out','ylim',[-1100 100]/1e3)
set(gca,'fontsize',14)
ylabel('displacement of maximum reflectivity below cloud base')
xlabel('normalized frequency')
% print('-depsc',[put_gordon 'CMR_below_cloud_hist.eps']);

% CFAD of CMR(1 min) and height relative to cloud base (zht-z10)
rfledge=[-Inf -40:30 Inf];
zedge=[-Inf -1600:25:1000 Inf];
ct=zeros(length(rfledge),length(zedge));
for ir=1:length(rfledge)-1
    iir=zmx>=rfledge(ir) & zmx<rfledge(ir+1);
    if ~isempty(iir)
        ct(ir,:)=histc(zht(iir)-z10(iir),zedge);
    end
end
iir=zmx>=rfledge(end);
if ~isempty(iir)
    ct(length(rfledge),:)=histc(zht(iir)-z10(iir),zedge);
end

clf
imagesc(rfledge(2:end-1)+.5,(zedge(2:end-1)+12.5)/1e3,log2(ct(2:end-1,2:end-1)')); axis xy
hold on
plot(zmx,(zht-z10)/1e3,'r.','markersize',1)
contour(rfledge(2:end-1)+.5,(zedge(2:end-1)+12.5)/1e3,log2(ct(2:end-1,2:end-1)'),2:2:8,'k');
colormap(1-gray(10)); caxis([0 10])
hc=colorbar; set(hc,'ytick',0:10,'yticklabel',[0 2.^(1:10)])
set(get(hc,'ylabel'),'string','count','fontsize',14); 
axis([-33.5 20 -1.6 1.0])
set(gca,'fontsize',14)
xlabel('maximum reflectivity (dBZ)')
ylabel('height displacement from cloud base')
print('-depsc',[put_gordon 'joint_CMR_zdisp_hist.eps'])

% composite above and below CB-100m (35 percentile)
wedges=-5:.1:5;
n=zeros(3,length(wedges),length(w.height));
m=double(0*w.cfad(1:3,:,:));

% zht is height of CMR, z10 is ceilo cloud base
ii=zht-z10>-100 & zht-z10<0; % 50% of minutes find CMR 0-100 m below cloud base
m(1,:,:)=sum(double(w.cfad(ii,:,:)));
for zi=1:length(w.height)
    n(1,:,zi)=hist(w.mean(ii,zi),wedges);
end
ii=zht-z10<-100; % 35% CMR >100m below base
m(2,:,:)=sum(double(w.cfad(ii,:,:)));
for zi=1:length(w.height)
    n(2,:,zi)=hist(w.mean(ii,zi),wedges);
end
ii=zht-z10>0; % 15% CMR above cloud base
m(3,:,:)=sum(double(w.cfad(ii,:,:)));
for zi=1:length(w.height)
    n(3,:,zi)=hist(w.mean(ii,zi),wedges);
end


% % shade "clouds" contour "drizzle"
% subplot(2,1,1)
% contourf(wedges,w.height/1e3,squeeze(n(1,:,:))',6,'edgecolor','none')
% hold on
% contour(wedges,w.height/1e3,squeeze(n(2,:,:))',6,'k')
% axis([-2.5 2.5 .15 1.8])
% set(gca,'fontsize',14)
% xlabel('vertical velocity (m s^{-1})')
% ylabel('height (km)')

% composites of 1 min histograms
clf
contourf(w.bins(2:end-1)+0.1,w.height/1e3,squeeze(m(1,:,2:end-1)+m(3,:,2:end-1)),7,'edgecolor','none')
hold on
contour(w.bins(2:end-1)+0.1,w.height/1e3,squeeze(m(2,:,2:end-1)),6,'k')
axis([-2.5 2.5 .15 1.8])
set(gca,'fontsize',14)
xlabel('vertical velocity (m s^{-1})')
ylabel('height (km)')
%print('-depsc',[put_gordon 'cmposit_w_on_CMRdispl.eps'])

clf
subplot(2,3,1)
contourf(Z.bins(2:end-1)+1,Z.height/1e3,squeeze(sum(Z.cfad(zht-z10>-100,:,2:end-1))),1e5*(.5:.5:4),'edgecolor','none')
hold on
contour(Z.bins(2:end-1)+1,Z.height/1e3,squeeze(sum(Z.cfad(zht-z10<-100,:,2:end-1))),1e5*(.3:.3:1.8),'color',[0 0 .6])
axis([-33 0 .15 1.8])
caxis([0 4e5])
set(gca,'fontsize',14)
xlabel('reflectivity (dBZ)')
ylabel('height (km)')

subplot(2,3,2)
contourf(w.bins(2:end-1)+0.1,w.height/1e3,squeeze(sum(w.cfad(zht-z10>-100,:,2:end-1))),7,'edgecolor','none')
hold on
contour(w.bins(2:end-1)+0.1,w.height/1e3,squeeze(sum(w.cfad(zht-z10<-100,:,2:end-1))),6,'color',[0 0 .6])
plot([0 0],[0 2],'k:')
axis([-2.5 2.5 .15 1.8])
set(gca,'fontsize',14)
xlabel({'Doppler vertical velocity (m s^{-1})'})

subplot(2,3,3)
contourf(Dw.bins(2:end-1)+0.1,Dw.height/1e3,squeeze(sum(Dw.cfad(zht-z10>-100,:,2:end-1))),7,'edgecolor','none')
hold on
contour(Dw.bins(2:end-1)+0.1,Dw.height/1e3,squeeze(sum(Dw.cfad(zht-z10<-100,:,2:end-1))),6,'color',[0 0 .6])
axis([0 1.5 .15 1.8])
set(gca,'fontsize',14)
xlabel({'Doppler velocity width (m s^{-1})'})
orient landscape
%print('-depsc',[put_gordon 'cmposit_ZwDw_onCMRdispl.eps'])
