load('./retrieval/Pinsky_lookup_allVOCALS_cloudrelative.mat')
Vgbin=squeeze(nanmean(Vgbin,3));
if ~exist('h','var')
    h=(33.7170:24.9756:3006)';
end
hrel=h(1:101)-h(91);

clf
pcolor(Zbin,hrel,Vgbin'); shading flat
hold on
contour(Zbin+0.125,hrel+12.5,Vgbin',[0 0],'edgecolor',[0.7 0 0])
contour(Zbin+0.125,hrel+12.5,squeeze(sum(nbin,2))',4.^(1:8),'edgecolor',[0 0 0])
contour(Zbin+0.125,hrel+12.5,squeeze(sum(nbin,2))',[30 30],'edgecolor',[0 1 0])
plot([-40 30],[0 0],'k--')
set(gca,'fontsize',14)
ylabel('height-cloud top height (m)')
set(gca,'ylim',hrel([20 end]),'xlim',[-40 30],'color',0.7+[0 0 0])
xlabel('dBZ bin')
colorbar
caxis([-1 1])
title('Vg(h,Z) lookup table')

%% difference plot
pcolor(Zbin,hrel,Vgbin'-Vgbin_Z_h'); shading flat
