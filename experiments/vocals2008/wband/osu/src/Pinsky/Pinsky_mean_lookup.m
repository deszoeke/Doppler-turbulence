% deprecated

% concatenate all the lookup tables
% cd /Users/sdeszoek/Data/cruises/VOCALS_2008/RHB/Scientific_analysis/programs/VOCALS2008_programs/wband/Pinsky/retrieval
% dr=dir('./');
% sprintf('''%s'',',dr(3:end).name)
% MATrcat(3,[3 4],'Pinsky_lookup2008_318_15.mat','Pinsky_lookup2008_318_16.mat','Pinsky_lookup2008_318_17.mat','Pinsky_lookup2008_318_18.mat','Pinsky_lookup2008_318_19.mat','Pinsky_lookup2008_318_20.mat','Pinsky_lookup2008_318_21.mat','Pinsky_lookup2008_318_22.mat','Pinsky_lookup2008_318_23.mat','Pinsky_lookup2008_319_00.mat','Pinsky_lookup2008_319_01.mat','Pinsky_lookup2008_319_02.mat','Pinsky_lookup2008_319_03.mat','Pinsky_lookup2008_319_04.mat','Pinsky_lookup2008_319_05.mat','Pinsky_lookup2008_319_06.mat','Pinsky_lookup2008_319_07.mat','Pinsky_lookup2008_319_08.mat','Pinsky_lookup2008_319_09.mat','Pinsky_lookup2008_319_10.mat','Pinsky_lookup2008_319_11.mat','Pinsky_lookup2008_319_12.mat','Pinsky_lookup2008_319_13.mat','Pinsky_lookup2008_319_14.mat','Pinsky_lookup2008_319_15.mat','Pinsky_lookup2008_319_16.mat','Pinsky_lookup2008_319_17.mat','Pinsky_lookup2008_319_18.mat','Pinsky_lookup2008_319_19.mat','Pinsky_lookup2008_319_20.mat','Pinsky_lookup2008_319_21.mat','Pinsky_lookup2008_319_22.mat','Pinsky_lookup2008_319_23.mat','Pinsky_lookup2008_320_00.mat','Pinsky_lookup2008_320_01.mat','Pinsky_lookup2008_320_02.mat','Pinsky_lookup2008_320_03.mat','Pinsky_lookup2008_320_04.mat','Pinsky_lookup2008_320_05.mat','Pinsky_lookup2008_320_06.mat','Pinsky_lookup2008_320_07.mat','Pinsky_lookup2008_320_08.mat','Pinsky_lookup2008_320_09.mat','Pinsky_lookup2008_320_10.mat','Pinsky_lookup2008_320_11.mat','Pinsky_lookup2008_320_12.mat','Pinsky_lookup2008_320_13.mat','Pinsky_lookup2008_320_14.mat','Pinsky_lookup2008_320_15.mat','Pinsky_lookup2008_320_16.mat','Pinsky_lookup2008_320_17.mat','Pinsky_lookup2008_320_18.mat','Pinsky_lookup2008_320_19.mat','Pinsky_lookup2008_320_20.mat','Pinsky_lookup2008_320_21.mat','Pinsky_lookup2008_320_22.mat','Pinsky_lookup2008_320_23.mat','Pinsky_lookup2008_321_00.mat','Pinsky_lookup2008_321_01.mat','Pinsky_lookup2008_321_02.mat','Pinsky_lookup2008_321_03.mat','Pinsky_lookup2008_321_04.mat','Pinsky_lookup2008_321_05.mat','Pinsky_lookup2008_321_06.mat','Pinsky_lookup2008_321_07.mat','Pinsky_lookup2008_321_08.mat','Pinsky_lookup2008_321_09.mat','Pinsky_lookup2008_321_10.mat','Pinsky_lookup2008_321_11.mat','Pinsky_lookup2008_321_12.mat','Pinsky_lookup2008_321_13.mat','Pinsky_lookup2008_321_14.mat','Pinsky_lookup2008_321_15.mat','Pinsky_lookup2008_321_16.mat','Pinsky_lookup2008_321_17.mat','Pinsky_lookup2008_321_18.mat','Pinsky_lookup2008_321_19.mat','Pinsky_lookup2008_321_20.mat','Pinsky_lookup2008_321_21.mat','Pinsky_lookup2008_321_22.mat','Pinsky_lookup2008_321_23.mat','Pinsky_lookup2008_322_00.mat','Pinsky_lookup2008_322_01.mat','Pinsky_lookup2008_322_02.mat','Pinsky_lookup2008_322_03.mat','Pinsky_lookup2008_322_04.mat','Pinsky_lookup2008_322_05.mat','Pinsky_lookup2008_322_06.mat','Pinsky_lookup2008_322_07.mat','Pinsky_lookup2008_322_08.mat','Pinsky_lookup2008_322_09.mat','Pinsky_lookup2008_322_10.mat','Pinsky_lookup2008_322_11.mat','Pinsky_lookup2008_322_12.mat','Pinsky_lookup2008_322_13.mat','Pinsky_lookup2008_322_14.mat','Pinsky_lookup2008_322_15.mat','Pinsky_lookup2008_322_16.mat','Pinsky_lookup2008_322_17.mat','Pinsky_lookup2008_322_18.mat','Pinsky_lookup2008_322_19.mat','Pinsky_lookup2008_322_20.mat','Pinsky_lookup2008_322_21.mat','Pinsky_lookup2008_322_22.mat','Pinsky_lookup2008_322_23.mat','Pinsky_lookup2008_323_00.mat','Pinsky_lookup2008_323_01.mat','Pinsky_lookup2008_323_02.mat','Pinsky_lookup2008_323_03.mat','Pinsky_lookup2008_323_04.mat','Pinsky_lookup2008_323_05.mat','Pinsky_lookup2008_323_06.mat','Pinsky_lookup2008_323_07.mat','Pinsky_lookup2008_323_08.mat','Pinsky_lookup2008_323_09.mat','Pinsky_lookup2008_323_10.mat','Pinsky_lookup2008_323_11.mat','Pinsky_lookup2008_323_12.mat','Pinsky_lookup2008_323_13.mat','Pinsky_lookup2008_323_14.mat','Pinsky_lookup2008_323_15.mat','Pinsky_lookup2008_323_16.mat','Pinsky_lookup2008_323_17.mat','Pinsky_lookup2008_323_18.mat','Pinsky_lookup2008_323_19.mat','Pinsky_lookup2008_323_20.mat','Pinsky_lookup2008_323_21.mat','Pinsky_lookup2008_323_22.mat','Pinsky_lookup2008_323_23.mat','Pinsky_lookup2008_324_00.mat','Pinsky_lookup2008_324_01.mat','Pinsky_lookup2008_324_02.mat','Pinsky_lookup2008_324_03.mat','Pinsky_lookup2008_324_04.mat','Pinsky_lookup2008_324_05.mat','Pinsky_lookup2008_324_06.mat','Pinsky_lookup2008_324_07.mat','Pinsky_lookup2008_324_08.mat','Pinsky_lookup2008_324_09.mat','Pinsky_lookup2008_324_10.mat','Pinsky_lookup2008_324_11.mat','Pinsky_lookup2008_324_12.mat','Pinsky_lookup2008_324_13.mat','Pinsky_lookup2008_324_14.mat','Pinsky_lookup2008_324_15.mat','Pinsky_lookup2008_324_16.mat','Pinsky_lookup2008_324_17.mat','Pinsky_lookup2008_324_18.mat','Pinsky_lookup2008_324_19.mat','Pinsky_lookup2008_324_20.mat','Pinsky_lookup2008_324_21.mat','Pinsky_lookup2008_324_22.mat','Pinsky_lookup2008_324_23.mat','Pinsky_lookup2008_325_00.mat','Pinsky_lookup2008_325_01.mat','Pinsky_lookup2008_325_02.mat','Pinsky_lookup2008_325_03.mat','Pinsky_lookup2008_325_04.mat','Pinsky_lookup2008_325_05.mat','Pinsky_lookup2008_325_06.mat','Pinsky_lookup2008_325_07.mat','Pinsky_lookup2008_325_08.mat','Pinsky_lookup2008_325_09.mat','Pinsky_lookup2008_325_10.mat','Pinsky_lookup2008_325_11.mat','Pinsky_lookup2008_325_12.mat','Pinsky_lookup2008_325_13.mat','Pinsky_lookup2008_325_14.mat','Pinsky_lookup2008_325_15.mat','Pinsky_lookup2008_325_16.mat','Pinsky_lookup2008_325_17.mat','Pinsky_lookup2008_325_18.mat','Pinsky_lookup2008_325_19.mat','Pinsky_lookup2008_325_20.mat','Pinsky_lookup2008_325_21.mat','Pinsky_lookup2008_325_22.mat','Pinsky_lookup2008_325_23.mat','Pinsky_lookup2008_326_00.mat','Pinsky_lookup2008_326_01.mat','Pinsky_lookup2008_326_02.mat','Pinsky_lookup2008_326_03.mat','Pinsky_lookup2008_326_04.mat','Pinsky_lookup2008_326_05.mat','Pinsky_lookup2008_326_06.mat','Pinsky_lookup2008_326_07.mat','Pinsky_lookup2008_326_08.mat','Pinsky_lookup2008_326_09.mat','Pinsky_lookup2008_326_10.mat','Pinsky_lookup2008_326_11.mat','Pinsky_lookup2008_326_12.mat','Pinsky_lookup2008_326_13.mat','Pinsky_lookup2008_326_14.mat','Pinsky_lookup2008_326_15.mat','Pinsky_lookup2008_326_16.mat','Pinsky_lookup2008_326_17.mat','Pinsky_lookup2008_326_18.mat','Pinsky_lookup2008_326_19.mat','Pinsky_lookup2008_326_20.mat','Pinsky_lookup2008_326_21.mat','Pinsky_lookup2008_326_22.mat','Pinsky_lookup2008_326_23.mat','Pinsky_lookup2008_327_00.mat','Pinsky_lookup2008_327_01.mat','Pinsky_lookup2008_327_02.mat','Pinsky_lookup2008_327_03.mat','Pinsky_lookup2008_327_04.mat','Pinsky_lookup2008_327_05.mat','Pinsky_lookup2008_327_06.mat','Pinsky_lookup2008_327_07.mat','Pinsky_lookup2008_327_08.mat','Pinsky_lookup2008_327_09.mat','Pinsky_lookup2008_327_10.mat','Pinsky_lookup2008_327_11.mat','Pinsky_lookup2008_327_12.mat','Pinsky_lookup2008_327_13.mat','Pinsky_lookup2008_327_14.mat','Pinsky_lookup2008_327_15.mat','Pinsky_lookup2008_327_16.mat','Pinsky_lookup2008_327_17.mat','Pinsky_lookup2008_327_18.mat','Pinsky_lookup2008_327_19.mat','Pinsky_lookup2008_327_20.mat','Pinsky_lookup2008_327_21.mat','Pinsky_lookup2008_327_22.mat','Pinsky_lookup2008_327_23.mat','Pinsky_lookup2008_328_00.mat','Pinsky_lookup2008_328_01.mat','Pinsky_lookup2008_328_02.mat','Pinsky_lookup2008_328_03.mat','Pinsky_lookup2008_328_04.mat','Pinsky_lookup2008_328_05.mat','Pinsky_lookup2008_328_06.mat','Pinsky_lookup2008_328_07.mat','Pinsky_lookup2008_328_08.mat','Pinsky_lookup2008_328_09.mat','Pinsky_lookup2008_328_10.mat','Pinsky_lookup2008_328_11.mat','Pinsky_lookup2008_328_12.mat','Pinsky_lookup2008_328_13.mat','Pinsky_lookup2008_328_14.mat','Pinsky_lookup2008_328_15.mat','Pinsky_lookup2008_328_16.mat','Pinsky_lookup2008_328_17.mat','Pinsky_lookup2008_328_18.mat','Pinsky_lookup2008_328_19.mat','Pinsky_lookup2008_328_20.mat','Pinsky_lookup2008_328_21.mat','Pinsky_lookup2008_328_22.mat','Pinsky_lookup2008_328_23.mat','Pinsky_lookup2008_329_00.mat','Pinsky_lookup2008_329_01.mat','Pinsky_lookup2008_329_02.mat','Pinsky_lookup2008_329_03.mat','Pinsky_lookup2008_329_04.mat','Pinsky_lookup2008_329_05.mat','Pinsky_lookup2008_329_06.mat','Pinsky_lookup2008_329_07.mat','Pinsky_lookup2008_329_08.mat','Pinsky_lookup2008_329_09.mat','Pinsky_lookup2008_329_10.mat','Pinsky_lookup2008_329_11.mat','Pinsky_lookup2008_329_12.mat','Pinsky_lookup2008_329_13.mat','Pinsky_lookup2008_329_14.mat','Pinsky_lookup2008_329_15.mat','Pinsky_lookup2008_329_16.mat','Pinsky_lookup2008_329_17.mat','Pinsky_lookup2008_329_18.mat','Pinsky_lookup2008_329_19.mat','Pinsky_lookup2008_329_20.mat','Pinsky_lookup2008_329_21.mat','Pinsky_lookup2008_329_22.mat','Pinsky_lookup2008_329_23.mat','Pinsky_lookup2008_330_00.mat','Pinsky_lookup2008_330_01.mat','Pinsky_lookup2008_330_02.mat','Pinsky_lookup2008_330_03.mat','Pinsky_lookup2008_330_04.mat','Pinsky_lookup2008_330_05.mat','Pinsky_lookup2008_330_06.mat','Pinsky_lookup2008_330_07.mat','Pinsky_lookup2008_330_08.mat','Pinsky_lookup2008_330_09.mat','Pinsky_lookup2008_330_10.mat','Pinsky_lookup2008_330_11.mat','Pinsky_lookup2008_330_12.mat','Pinsky_lookup2008_330_13.mat','Pinsky_lookup2008_330_14.mat','Pinsky_lookup2008_330_15.mat','Pinsky_lookup2008_330_16.mat','Pinsky_lookup2008_330_17.mat','Pinsky_lookup2008_330_18.mat','Pinsky_lookup2008_330_19.mat','Pinsky_lookup2008_330_20.mat','Pinsky_lookup2008_330_21.mat','Pinsky_lookup2008_330_22.mat','Pinsky_lookup2008_330_23.mat','Pinsky_lookup2008_331_00.mat','Pinsky_lookup2008_331_01.mat','Pinsky_lookup2008_331_02.mat','Pinsky_lookup2008_331_03.mat','Pinsky_lookup2008_331_04.mat','Pinsky_lookup2008_331_05.mat','Pinsky_lookup2008_331_06.mat','Pinsky_lookup2008_331_07.mat','Pinsky_lookup2008_331_08.mat','Pinsky_lookup2008_331_09.mat','Pinsky_lookup2008_331_10.mat','Pinsky_lookup2008_331_11.mat','Pinsky_lookup2008_331_12.mat','Pinsky_lookup2008_331_13.mat','Pinsky_lookup2008_331_14.mat','Pinsky_lookup2008_331_15.mat','Pinsky_lookup2008_331_16.mat','Pinsky_lookup2008_331_17.mat','Pinsky_lookup2008_331_18.mat','Pinsky_lookup2008_331_19.mat','Pinsky_lookup2008_331_20.mat','Pinsky_lookup2008_331_21.mat','Pinsky_lookup2008_331_22.mat','Pinsky_lookup2008_331_23.mat','Pinsky_lookup2008_332_00.mat','Pinsky_lookup2008_332_01.mat','Pinsky_lookup2008_332_02.mat','Pinsky_lookup2008_332_03.mat','Pinsky_lookup2008_332_04.mat','Pinsky_lookup2008_332_05.mat','Pinsky_lookup2008_332_06.mat','Pinsky_lookup2008_332_07.mat','Pinsky_lookup2008_332_08.mat','Pinsky_lookup2008_332_09.mat','Pinsky_lookup2008_332_10.mat','Pinsky_lookup2008_332_11.mat','Pinsky_lookup2008_332_12.mat','Pinsky_lookup2008_332_13.mat','Pinsky_lookup2008_332_14.mat','Pinsky_lookup2008_332_15.mat','Pinsky_lookup2008_332_16.mat','Pinsky_lookup2008_332_17.mat','Pinsky_lookup2008_332_18.mat','Pinsky_lookup2008_332_19.mat','Pinsky_lookup2008_332_20.mat','Pinsky_lookup2008_332_21.mat','Pinsky_lookup2008_332_22.mat','Pinsky_lookup2008_332_23.mat','Pinsky_lookup2008_333_00.mat','Pinsky_lookup2008_333_01.mat','Pinsky_lookup2008_333_02.mat','Pinsky_lookup2008_333_03.mat','Pinsky_lookup2008_333_04.mat','Pinsky_lookup2008_333_05.mat','Pinsky_lookup2008_333_06.mat','Pinsky_lookup2008_333_07.mat','Pinsky_lookup2008_333_08.mat','Pinsky_lookup2008_333_09.mat','Pinsky_lookup2008_333_10.mat','Pinsky_lookup2008_333_11.mat','Pinsky_lookup2008_333_12.mat','Pinsky_lookup2008_333_13.mat','Pinsky_lookup2008_333_14.mat','Pinsky_lookup2008_333_15.mat','Pinsky_lookup2008_333_16.mat','Pinsky_lookup2008_333_17.mat','Pinsky_lookup2008_333_18.mat','Pinsky_lookup2008_333_19.mat','Pinsky_lookup2008_333_20.mat','Pinsky_lookup2008_333_21.mat','Pinsky_lookup2008_333_22.mat','Pinsky_lookup2008_333_23.mat','Pinsky_lookup2008_334_00.mat','Pinsky_lookup2008_334_01.mat','Pinsky_lookup2008_334_02.mat','Pinsky_lookup2008_334_03.mat','Pinsky_lookup2008_334_04.mat','Pinsky_lookup2008_334_05.mat','Pinsky_lookup2008_334_06.mat','Pinsky_lookup2008_334_07.mat','Pinsky_lookup2008_334_08.mat','Pinsky_lookup2008_334_09.mat','Pinsky_lookup2008_334_10.mat','Pinsky_lookup2008_334_11.mat','Pinsky_lookup2008_334_12.mat','Pinsky_lookup2008_334_13.mat','Pinsky_lookup2008_334_14.mat','Pinsky_lookup2008_334_15.mat','Pinsky_lookup2008_334_16.mat','Pinsky_lookup2008_334_17.mat','Pinsky_lookup2008_334_18.mat','Pinsky_lookup2008_334_19.mat','Pinsky_lookup2008_334_20.mat','Pinsky_lookup2008_334_21.mat','Pinsky_lookup2008_334_22.mat','Pinsky_lookup2008_334_23.mat','Pinsky_lookup2008_335_00.mat','Pinsky_lookup2008_335_01.mat','Pinsky_lookup2008_335_02.mat','Pinsky_lookup2008_335_03.mat','Pinsky_lookup2008_335_04.mat','Pinsky_lookup2008_335_05.mat','Pinsky_lookup2008_335_06.mat','Pinsky_lookup2008_335_07.mat','Pinsky_lookup2008_335_08.mat','Pinsky_lookup2008_335_09.mat','Pinsky_lookup2008_335_10.mat','Pinsky_lookup2008_335_11.mat','Pinsky_lookup2008_335_12.mat','Pinsky_lookup2008_335_13.mat','Pinsky_lookup2008_335_14.mat','Pinsky_lookup2008_335_15.mat','Pinsky_lookup2008_335_16.mat','Pinsky_lookup2008_335_17.mat','Pinsky_lookup2008_335_18.mat','Pinsky_lookup2008_335_19.mat','Pinsky_lookup2008_335_20.mat','Pinsky_lookup2008_335_21.mat','Pinsky_lookup2008_335_22.mat','Pinsky_lookup2008_335_23.mat','Pinsky_lookup_allVOCALS.mat')

cd /Users/sdeszoek/Data/cruises/VOCALS_2008/RHB/Scientific_analysis/programs/VOCALS2008_programs/wband/Pinsky
clear
load('./retrieval/Pinsky_lookup_allVOCALS.mat')

% by eyeballing, retain 9 EOFs
neof=9;

countt=sum(isfinite(Vgbin),3); % count number of times for each h-Z point

b2rcolormap(17); dock
clf
subplot(2,1,1)
contourf(Zbin,h,nanmean(Vgbin,3)',-3:.2:.1,'edgecolor','none');
caxis([-3 0.2])
set(gca,'color',0.7*[1 1 1],'ylim',[150 2000],'fontsize',12)
colorbar
hold on
contour(Zbin,h,log2(countt'),1:9,'k');
contour(Zbin,h,log2(countt'),[3 3],'b');
title('mean <V>(h,Z)')

subplot(2,1,2)
contourf(Zbin,h,nanstd(Vgbin,3)',0:0.1:0.8,'edgecolor','none');
caxis([0 0.8])
colorbar
hold on
contour(Zbin,h,log2(countt'),1:9,'k');
contour(Zbin,h,log2(countt'),[3 3],'b');
set(gca,'color',0.7*[1 1 1],'ylim',[150 2000],'fontsize',12)
title('standard deviation among hours of <V>(h,Z)')
xlabel('dBZ')
ylabel('height (m)')

print -dpng Pinsky_mean_lookup.png

mesh=meshgrid(h',Zbin');
indx=countt>100; % collapse to parts of the h-Z space with statistically significantly sampled variability

A=reshape(Vgbin,[120*201,417]);
A=A(indx,:);
A=A(:,sum(isfinite(A))>100);
Ap=A-repmat(nanmean(A,2),[1 size(A,2)]);
np=sum(indx(:)); % structure dimension
nq=size(A,2);    % sampling (time) dimension

% compute the EOFs better handling sparseness by computing the covariance matrix ignoring NaNs
Cov=nancov(Ap');
[pe,se]=eigs(Cov,neof); % modes are h-Z structure vectors
sed=diag(se);
stdpe=std(pe); % for normalizing pe to have unit variance
Apf=Ap;
Apf(isnan(Apf))=0;
% dimensional time-varying amplitude of the modes
qe=Apf'*pe./(repmat(sum(isfinite(Ap))'-1,[1 neof]).*repmat(stdpe,[nq 1]));
% normalized nondimensional amplitude of the modes
qen=qe./std(qe);

%%%
stdAq=nanstd(Ap);  % structural standard deviation as a function of time
stdAp=nanstd(Ap,2); % time standard deviation for each point in h-Z space
qprod=Apf'*pe; % non normalized
qcov=qprod./repmat(sum(isfinite(Ap))'-1,[1 neof]);
qcor=qcov./(repmat(stdAq',[1 neof]).*repmat(std(pe),[nq 1]));

plot(sed,'.-')
Pe=NaN(size(indx,1),size(indx,2),neof);
[ii,jj]=find(indx);
for kk=1:neof;
    for ss=1:length(ii)
        Pe(ii(ss),jj(ss),kk)=P(ss,kk);
    end
end
for kk=1:neof
    subplot(3,3,kk)
    pcolor(Zbin(1:181),h(5:75)/1e3,sed(kk)*Pe(1:181,5:75,kk)'); shading flat
    colorbar('east')
    caxis(max(abs(caxis))*[-1 1])
end

% old SVD way unfortunately fills missing data with zeros
Ap(isnan(Ap))=0;
[P,S,Q]=svd(Ap);
% P is h&Z structure, Q is hourly sampling structure
s=diag(S);
% refill Z-h space
p=NaN(size(indx,1),size(indx,2),neof);
[ii,jj]=find(indx);
for kk=1:neof;
    for ss=1:length(ii)
        p(ii(ss),jj(ss),kk)=P(ss,kk);
    end
end
q=Q(1:neof,:);
for kk=1:neof
    subplot(3,3,kk)
    pcolor(Zbin(1:181),h(5:75)/1e3,s(kk)*p(1:181,5:75,kk)'); shading flat
    colorbar('east')
    caxis(max(abs(caxis))*[-1 1])
end