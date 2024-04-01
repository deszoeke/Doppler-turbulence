% call after Pinsky_retrieval or Pinsky_retrieval_cloudrelative

b2rcolormap(25);

if 0
    %% test plot the height dependent lookup
    clf
    subplot(2,1,1)
    pcolor(Zbin-0.125,h-15,Vgbin'); shading flat
    hold on
    contour(Zbin,h,Vgbin',0:0.1:1,'edgecolor',0.7+[0 0 0])
    set(gca,'fontsize',14)
    ylabel('height (m)')
    set(gca,'ylim',[0 1800],'color',0.7+[0 0 0])
    xlabel('dBZ bin')
    colorbar
    caxis([-1.5 0])
    title('Vg(h,Z) lookup table \phi')
    axis([-40 0 0 1600])
    
    subplot(2,1,2)
    pcolor(Zbin-0.125,h-15,Ustdbin'); shading flat
    set(gca,'fontsize',14)
    ylabel('height (m)')
    set(gca,'ylim',[0 1800],'color',0.7+[0 0 0])
    xlabel('dBZ bin')
    colorbar
    title('<U^2>^{1/2}(h,Z) lookup table \theta')
    axis([-40 0 0 1600])
    caxis([.2 1.2])
    orient tall
    % print -dpng Pinsky_lookup_tables.png
end

%% test plot the cloud-relative lookup
hrel=h(1:101)-h(91);
clf
subplot(2,1,1)
pcolor(Zbin-0.125,hrel-15,Vgbin'); shading flat
hold on
contour(Zbin,hrel,Vgbin',0:0.1:1,'edgecolor',0.7+[0 0 0])
set(gca,'fontsize',14)
ylabel('height-cloud top height (m)')
set(gca,'ylim',hrel([26 end-1]),'xlim',[-40 30],'color',0.7+[0 0 0])
xlabel('dBZ bin')
colorbar
caxis([-3 0])
title('Vg(h,Z) lookup table \phi')

subplot(2,1,2)
pcolor(Zbin-0.125,hrel-15,Ustdbin'); shading flat
set(gca,'fontsize',14)
ylabel('height-cloud top height (m)')
set(gca,'ylim',hrel([26 end]),'xlim',[-40 30],'color',0.7+[0 0 0])
xlabel('dBZ bin')
colorbar
title('<U^2>^{1/2}(h,Z) lookup table \theta')
caxis([.1 1.2])

%% set up cloud-relative lookup for all hours/days
hrel=h(1:101)-h(91);
cd retrieval 
MATrcat(3,[3 4 5],'Pinsky_lookup2008_318_15_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_318_16_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_318_17_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_318_18_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_318_19_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_318_20_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_318_21_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_318_22_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_318_23_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_319_00_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_319_01_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_319_02_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_319_03_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_319_04_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_319_05_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_319_06_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_319_07_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_319_08_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_319_09_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_319_10_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_319_11_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_319_12_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_319_13_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_319_14_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_319_15_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_319_16_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_319_17_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_319_18_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_319_19_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_319_20_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_319_21_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_319_22_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_319_23_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_320_00_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_320_01_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_320_02_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_320_03_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_320_04_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_320_05_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_320_06_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_320_07_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_320_08_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_320_09_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_320_10_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_320_11_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_320_12_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_320_13_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_320_14_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_320_15_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_320_16_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_320_17_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_320_18_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_320_19_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_320_20_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_320_21_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_320_22_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_320_23_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_321_00_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_321_01_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_321_02_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_321_03_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_321_04_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_321_05_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_321_06_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_321_07_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_321_08_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_321_09_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_321_10_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_321_11_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_321_12_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_321_13_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_321_14_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_321_15_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_321_16_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_321_17_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_321_18_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_321_19_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_321_20_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_321_21_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_321_22_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_321_23_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_322_00_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_322_01_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_322_02_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_322_03_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_322_04_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_322_05_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_322_06_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_322_07_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_322_08_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_322_09_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_322_10_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_322_11_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_322_12_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_322_13_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_322_14_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_322_15_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_322_16_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_322_17_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_322_18_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_322_19_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_322_20_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_322_21_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_322_22_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_322_23_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_323_00_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_323_01_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_323_02_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_323_03_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_323_04_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_323_05_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_323_06_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_323_07_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_323_08_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_323_09_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_323_10_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_323_11_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_323_12_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_323_13_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_323_14_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_323_15_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_323_16_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_323_17_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_323_18_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_323_19_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_323_20_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_323_21_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_323_22_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_323_23_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_324_00_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_324_01_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_324_02_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_324_03_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_324_04_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_324_05_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_324_06_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_324_07_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_324_08_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_324_09_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_324_10_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_324_11_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_324_12_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_324_13_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_324_14_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_324_15_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_324_16_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_324_17_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_324_18_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_324_19_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_324_20_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_324_21_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_324_22_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_324_23_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_325_00_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_325_01_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_325_02_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_325_03_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_325_04_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_325_05_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_325_06_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_325_07_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_325_08_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_325_09_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_325_10_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_325_11_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_325_12_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_325_13_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_325_14_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_325_15_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_325_16_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_325_17_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_325_18_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_325_19_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_325_20_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_325_21_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_325_22_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_325_23_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_326_00_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_326_01_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_326_02_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_326_03_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_326_04_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_326_05_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_326_06_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_326_07_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_326_08_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_326_09_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_326_10_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_326_11_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_326_12_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_326_13_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_326_14_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_326_15_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_326_16_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_326_17_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_326_18_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_326_19_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_326_20_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_326_21_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_326_22_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_326_23_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_327_00_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_327_01_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_327_02_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_327_03_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_327_04_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_327_05_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_327_06_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_327_07_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_327_08_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_327_09_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_327_10_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_327_11_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_327_12_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_327_13_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_327_14_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_327_15_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_327_16_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_327_17_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_327_18_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_327_19_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_327_20_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_327_21_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_327_22_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_327_23_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_328_00_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_328_01_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_328_02_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_328_03_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_328_04_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_328_05_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_328_06_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_328_07_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_328_08_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_328_09_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_328_10_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_328_11_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_328_12_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_328_13_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_328_14_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_328_15_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_328_16_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_328_17_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_328_18_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_328_19_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_328_20_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_328_21_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_328_22_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_328_23_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_329_00_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_329_01_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_329_02_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_329_03_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_329_04_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_329_05_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_329_06_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_329_07_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_329_08_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_329_09_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_329_10_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_329_11_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_329_12_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_329_13_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_329_14_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_329_15_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_329_16_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_329_17_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_329_18_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_329_19_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_329_20_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_329_21_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_329_22_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_329_23_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_330_00_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_330_01_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_330_02_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_330_03_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_330_04_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_330_05_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_330_06_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_330_07_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_330_08_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_330_09_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_330_10_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_330_11_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_330_12_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_330_13_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_330_14_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_330_15_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_330_16_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_330_17_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_330_18_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_330_19_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_330_20_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_330_21_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_330_22_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_330_23_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_331_00_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_331_01_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_331_02_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_331_03_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_331_04_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_331_05_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_331_06_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_331_07_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_331_08_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_331_09_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_331_10_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_331_11_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_331_12_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_331_13_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_331_14_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_331_15_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_331_16_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_331_17_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_331_18_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_331_19_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_331_20_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_331_21_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_331_22_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_331_23_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_332_00_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_332_01_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_332_02_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_332_03_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_332_04_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_332_05_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_332_06_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_332_07_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_332_08_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_332_09_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_332_10_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_332_11_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_332_12_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_332_13_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_332_14_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_332_15_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_332_16_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_332_17_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_332_18_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_332_19_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_332_20_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_332_21_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_332_22_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_332_23_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_333_00_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_333_01_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_333_02_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_333_03_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_333_04_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_333_05_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_333_06_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_333_07_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_333_08_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_333_09_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_333_10_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_333_11_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_333_12_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_333_13_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_333_14_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_333_15_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_333_16_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_333_17_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_333_18_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_333_19_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_333_20_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_333_21_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_333_22_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_333_23_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_334_00_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_334_01_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_334_02_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_334_03_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_334_04_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_334_05_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_334_06_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_334_07_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_334_08_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_334_09_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_334_10_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_334_11_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_334_12_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_334_13_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_334_14_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_334_15_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_334_16_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_334_17_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_334_18_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_334_19_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_334_20_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_334_21_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_334_22_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_334_23_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_335_00_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_335_01_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_335_02_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_335_03_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_335_04_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_335_05_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_335_06_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_335_07_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_335_08_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_335_09_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_335_10_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_335_11_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_335_12_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_335_13_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_335_14_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_335_15_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_335_16_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_335_17_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_335_18_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_335_19_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_335_20_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_335_21_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_335_22_cloudrelative_p5dBZ.mat','Pinsky_lookup2008_335_23_cloudrelative_p5dBZ.mat','Pinsky_lookup_allVOCALS_cloudrelative_p5dBZ.mat')
cd ..
load ./retrieval/Pinsky_lookup_allVOCALS_cloudrelative_p5dBZ.mat

%% plot all-days composite cloud relative lookup
clf
subplot(2,1,1)
pcolor(Zbin-0.25,hrel-15,nanmean(Vgbin,3)'); shading flat
hold on
contour(Zbin,hrel,nanmean(Vgbin,3)',[0 0],'edgecolor',[0.7 0 0])
contour(Zbin,hrel,sum(isfinite(Vgbin),3)',2.^(1:8),'edgecolor',[0 0 0])
contour(Zbin,hrel,sum(isfinite(Vgbin),3)',[8 8],'edgecolor',[0 1 0])
plot([-40 30],[0 0],'k--')
set(gca,'fontsize',14)
ylabel('height-cloud top height (m)')
set(gca,'ylim',hrel([20 end]),'xlim',[-40 30],'color',0.7+[0 0 0])
xlabel('dBZ bin')
colorbar
caxis([-3 1.5])
title('Vg(h,Z) lookup table \phi')

subplot(2,1,2)
pcolor(Zbin-0.25,hrel-15,nanmean(Ustdbin,3)'); shading flat
hold on
contour(Zbin,hrel,sum(isfinite(Vgbin),3)',2.^(1:8),'edgecolor',[0 0 0])
contour(Zbin,hrel,sum(isfinite(Vgbin),3)',[8 8],'edgecolor',[0 1 0])
plot([-40 30],[0 0],'k--')
set(gca,'fontsize',14)
ylabel('height-cloud top height (m)')
set(gca,'ylim',hrel([20 end]),'xlim',[-40 30],'color',0.7+[0 0 0])
xlabel('dBZ bin')
colorbar
title('<U^2>^{1/2}(h,Z) lookup table \theta')
caxis([.1 1.2])

orient tall
% print -dpng Pinsky_lookup_tables_all_cloudrelative.png

%% loop all retrieved hours
set(gcf,'renderer','zbuffer')
year=2008;
for iday=stday:enday;
    for idhr=00:23;
        fname=sprintf('/Users/sdeszoek/Data/cruises/VOCALS_2008/RHB/Scientific_analysis/programs/VOCALS2008_programs/wband/Pinsky/retrieval/Pinsky_wRetrieval%04i_%03i_%02i.mat',2008,iday,idhr);
        if exist(fname,'file')
            % load 1 hour
            load(fname);
            U=V-Vg;
            Vgprime=U-W;
            a=W./U;
            countt=length(time_offset);
            
            % plot all time-height series from 1 hour retrieval
            clf
            orient tall
            height=0.12;
            height_offset=0.015;
            width=0.775;
            left=0.125;
            bottom=0.87;
            ax(1)=axes('position',[left,bottom,width,height]);
            pcolor(time_offset(1:countt)/60,h,Z'); shading flat
            colorbar
            jtext('Z',0.02,0.16);
            ylabel('height (m)')
            set(gca,'ylim',[0 1600],'xtick',0:10:60,'xticklabel','');
            
            bottom=bottom-height-height_offset;
            ax(2)=axes('position',[left,bottom,width,height]);
            pcolor(time_offset(1:countt)/60,h,V'); shading flat
            colorbar
            caxis([-3 3])
            jtext('V',0.02,0.16);
            ylabel('height (m)')
            set(gca,'ylim',[0 1600],'xtick',0:10:60,'xticklabel','');
            
            bottom=bottom-height-height_offset;
            ax(3)=axes('position',[left,bottom,width,height]);
            pcolor(time_offset(1:countt)/60,h,Vg'); shading flat
            colorbar
            caxis([-2 2])
            jtext('Vg(h,Z(h,t))',0.02,0.16);
            ylabel('height (m)')
            set(gca,'ylim',[0 1600],'xtick',0:10:60,'xticklabel','');
            
            bottom=bottom-height-height_offset;
            ax(4)=axes('position',[left,bottom,width,height]);
            pcolor(time_offset(1:countt)/60,h,Ustd'); shading flat
            colorbar
            caxis([0.2 1.2])
            % caxis([0 0.7])
            jtext('\sigma_U(h,Z(h,t))',0.02,0.16);
            ylabel('height (m)')
            set(gca,'ylim',[0 1600],'xtick',0:10:60,'xticklabel','');
            
            bottom=bottom-height-height_offset;
            ax(5)=axes('position',[left,bottom,width,height]);
            pcolor(time_offset(1:countt)/60,h,U'); shading flat
            colorbar
            caxis([-1.5 1.5])
            jtext('residual velocity U',0.02,0.16);
            ylabel('height (m)')
            set(gca,'ylim',[0 1600],'xtick',0:10:60,'xticklabel','');
            
            bottom=bottom-height-height_offset;
            ax(6)=axes('position',[left,bottom,width,height]);
            pcolor(time_offset(1:countt)/60,h,Vgprime'); shading flat
            colorbar
            caxis([-.5 .5])
            jtext('Vg''',0.02,0.16);
            ylabel('height (m)')
            set(gca,'ylim',[0 1600],'xtick',0:10:60,'xticklabel','');
            
            bottom=bottom-height-height_offset;
            ax(7)=axes('position',[left,bottom,width,height]);
            pcolor(time_offset(1:countt)/60,h,W'); shading flat
            colorbar
            caxis([-1.5 1.5])
            jtext('W',0.02,0.16);
            ylabel('height (m)')
            set(gca,'ylim',[0 1600],'xtick',0:10:60);
            xlabel({'minute' sprintf('VOCALS %i day %i hour %i',year,iday,idhr)})
            set(ax(:),'color',0.7+[0 0 0])
            
            prname=sprintf('/Users/sdeszoek/Data/cruises/VOCALS_2008/RHB/Scientific_analysis/programs/VOCALS2008_programs/wband/Pinsky/retrieval/plot/Pinsky_wRetrieval_%03i_%02i.png',iday,idhr);
            print('-dpng',prname)
        end
    end
end

%% plot to 1-a0/sigma
axes(ax(4))
caxis([0 0.7])

axes(ax(5))
pcolor(time_offset(1:countt)/60,h,1-a'); shading flat
colorbar
caxis([-1 1])
jtext('1-a_0/\sigma_U',0.02,0.16);
ylabel('height (m)')
set(gca,'ylim',[0 1600],'xtick',0:10:60,'xticklabel','');

%% other plots
if 0
%% test plot the lookup
subplot(2,1,1)
plot(Zbin+0.125,Vgbin)
hold on
plot(Zbin+0.125,Ustdbin,'r')
plot(Zbin+0.125,log10(nbin)/10,'k:')
axis tight
set(gca,'ylim',[-1 1],'xlim',[-40 0])

%% test plot V and Vg
subplot(2,1,1)
imagesc(0:0.03:3600,h,V'); set(gca,'ydir','normal')
caxis([-3 3])
colorbar
vg=Vg;
vg(isnan(Vg))=999;
subplot(2,1,2)
imagesc(0:0.03:3600,h,vg'); set(gca,'ydir','normal')
colorbar
caxis([-1 0])
end