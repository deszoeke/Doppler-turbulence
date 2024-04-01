% run after Pinsky_main_cloudrelative.m
% cd /Users/sdeszoek/Data/cruises/VOCALS_2008/RHB/Scientific_analysis/programs/VOCALS2008_programs/wband/Pinsky

%% bins
Zbin=-40:0.25:30;
thickbin=0:100:1400;
if ~exist('h','var')
    h=(33.7170:24.9756:3006)';
end
hrel=h(1:101)-h(91);

%% record-average lookup tables

% composite whole-cruise 3D-bin lookup table from hourly tables
% Vg(dBZ,h-htop,cloudthickness)
stday=318;
enday=336;
[Vgbin,Ustdbin,nbin]=deal(zeros(length(Zbin),length(thickbin),length(hrel)));
for iday=stday:enday;
    for idhr=0:23;
        fname=sprintf('/Users/sdeszoek/Data/cruises/VOCALS_2008/RHB/Scientific_analysis/programs/VOCALS2008_programs/wband/Pinsky/retrieval/Pinsky_lookup%04i_%03i_%02i_cloudrelative_thickness.mat',2008,iday,idhr);
        if exist(fname,'file')
            Tmp=load(fname);
            Tmp.Vgbin(isnan(Tmp.Vgbin))=0;
            Tmp.Ustdbin(isnan(Tmp.Ustdbin))=0;
            Tmp.nbin(isnan(Tmp.Vgbin))=0;
            Vgbin=Vgbin+Tmp.Vgbin.*Tmp.nbin;
            Ustdbin=Ustdbin+Tmp.Ustdbin.*Tmp.nbin;
            nbin=nbin+Tmp.nbin;
        end
    end
end
Vgbin=Vgbin./nbin;
Ustdbin=Ustdbin./nbin;
Vgbin(nbin==0)=NaN;
Ustdbin(nbin==0)=NaN;

%% average over bin dimensions
% avg. thickness bins
Vgbin_Z_h=squeeze(nansum(Vgbin.*nbin,2)./sum(nbin,2));
Ustdbin_Z_h=squeeze(nansum(Ustdbin.*nbin,2)./sum(nbin,2));
% avg. dBZ bins
Vgbin_t_h=squeeze(nansum(Vgbin.*nbin,1)./sum(nbin,1));
Ustdbin_t_h=squeeze(nansum(Ustdbin.*nbin,1)./sum(nbin,1));
% avg. hrel bins
Vgbin_Z_t=squeeze(nansum(Vgbin.*nbin,3)./sum(nbin,3));
Ustdbin_Z_t=squeeze(nansum(Ustdbin.*nbin,3)./sum(nbin,3));

%% plot Vg(Z,h)
clf
pcolor(Zbin,hrel,Vgbin_Z_h'); shading flat
hold on
contour(Zbin+0.125,hrel+12.5,squeeze(Vgbin_Z_h)',[0 0],'edgecolor',[0.7 0 0])
contour(Zbin+0.125,hrel+12.5,squeeze(sum(nbin,2))',4.^(1:8),'edgecolor',[0 0 0])
contour(Zbin+0.125,hrel+12.5,squeeze(sum(nbin,2))',[30 30],'edgecolor',[0 1 0])
plot([-40 30],[0 0],'k--')
set(gca,'fontsize',14)
ylabel('height-cloud top height (m)')
set(gca,'ylim',hrel([20 end]),'xlim',[-40 30],'color',0.7+[0 0 0])
xlabel('dBZ bin')
colorbar
caxis([-4.5 1.5])
title('Vg(h,Z) lookup table')
%print -dpng Pinsky_lookup_Vg_Z_h.png

%% plot Vg(Z,t)
clf
pcolor(Zbin,thickbin,Vgbin_Z_t'); shading flat
hold on
contour(Zbin+0.125,thickbin+50,squeeze(Vgbin_Z_t)',[0 0],'edgecolor',[0.7 0 0])
contour(Zbin+0.125,thickbin+50,squeeze(sum(nbin,3))',4.^(1:8),'edgecolor',[0 0 0])
contour(Zbin+0.125,thickbin+50,squeeze(sum(nbin,3))',[30 30],'edgecolor',[0 1 0])
set(gca,'fontsize',14)
ylabel('cloud thickness (m)')
set(gca,'ylim',[0 1200],'xlim',[-40 25],'color',0.7+[0 0 0])
xlabel('dBZ bin')
colorbar
caxis([-3 1.5])
title('Vg(Z,thick) lookup table')
%print -dpng Pinsky_lookup_Vg_Z_t.png

%% plot Vg(t,h)
clf
pcolor(thickbin,hrel,Vgbin_t_h'); shading flat
hold on
contour(thickbin+50,hrel+12.5,squeeze(Vgbin_t_h)',[0 0],'edgecolor',[0.7 0 0])
contour(thickbin+50,hrel+12.5,squeeze(sum(nbin,1))',4.^(1:8),'edgecolor',[0 0 0])
contour(thickbin+50,hrel+12.5,squeeze(sum(nbin,1))',[30 30],'edgecolor',[0 1 0])
plot([-40 30],[0 0],'k--')
set(gca,'fontsize',14)
ylabel('height-cloud top height (m)')
set(gca,'ylim',hrel([20 end]),'xlim',[0 1200],'color',0.7+[0 0 0])
xlabel('cloud thickness (m)')
colorbar
caxis([-3 1.5])
title('Vg(thick,h) lookup table')
%print -dpng Pinsky_lookup_Vg_t_h.png

%% plot thickness bins separately
clf
offs=1;
for nth=1:6
    subplot(3,2,nth)
    pcolor(Zbin,hrel,squeeze(Vgbin(:,nth+offs,:))'); shading flat
    hold on
    contour(Zbin+0.125,hrel+12.5,squeeze(Vgbin(:,nth+offs,:))',[0 0],'edgecolor',[0.7 0 0])
    contour(Zbin+0.125,hrel+12.5,squeeze(nbin(:,nth+offs,:))',4.^(1:8),'edgecolor',[0 0 0])
    contour(Zbin+0.125,hrel+12.5,squeeze(nbin(:,nth+offs,:))',[30 30],'edgecolor',[0 1 0])
    plot([-40 30],[0 0],'k--')
    set(gca,'fontsize',14)
    ylabel('height-cloud top height (m)')
    set(gca,'ylim',hrel([20 end]),'xlim',[-40 30],'color',0.7+[0 0 0])
    xlabel('dBZ bin')
    colorbar
    caxis([-3 1.5])
    title(sprintf('Vg(h,Z) lookup table: cloud thickness %i-%i m', thickbin(nth+offs+(0:1))))
end
orient tall
%print -dpng Pinsky_lookup_Vg_all_cloudrelative_thickness.png

clf
for nth=1:6
    subplot(3,2,nth)
    pcolor(Zbin,hrel,squeeze(Ustdbin(:,nth+offs,:))'); shading flat
    hold on
    contour(Zbin+0.125,hrel+12.5,squeeze(nbin(:,nth+offs,:))',4.^(1:8),'edgecolor',[0 0 0])
    contour(Zbin+0.125,hrel+12.5,squeeze(nbin(:,nth+offs,:))',[30 30],'edgecolor',[0 1 0])
    plot([-40 30],[0 0],'k--')
    set(gca,'fontsize',14)
    ylabel('height-cloud top height (m)')
    set(gca,'ylim',hrel([20 end]),'xlim',[-40 30],'color',0.7+[0 0 0])
    xlabel('dBZ bin')
    colorbar
    title(sprintf('<U^2>^{1/2}(h,Z) lookup table: cloud thickness %i-%i m', thickbin(nth+offs+(0:1))));
    caxis([.1 1.2])
end
%print -dpng Pinsky_lookup_Ustd_all_cloudrelative_thickness.png

