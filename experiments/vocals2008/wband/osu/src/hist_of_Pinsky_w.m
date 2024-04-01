% hist_of_Pinsky_w.m
% Make 1-hour histogram of Pinsky air velocities

% cd /Volumes/Pegasus2/cruises/VOCALS_2008/RHB/Scientific_analysis/programs/VOCALS2008_programs/wband

datadir='/Volumes/Pegasus2/cruises/VOCALS_2008/RHB/Scientific_analysis/programs/VOCALS2008_programs/wband/Pinsky/retrieval';
filename=[datadir '/' 'Pinsky_wRetrieval2008_319_00.mat'];
load(filename);

% fractional area sampled
isf=isfinite(W);
nsampl=sum(isfinite(time_offset));
asampl=sum(isf,1)./nsampl;

wedges=(-7:0.5:7)';
count=histc(W,wedges); % carries an extra (empty) bin at end
acount=count./nsampl;
acondcount=bsxfun(@rdivide,count,asampl);

% mass flux
sw=sort(W);
csw=cumsum(sw);
[logi,ind]=max(sw>0); 
mid=zeros(size(ind));
for i=1:size(csw,2)
    mid(i)=csw(ind(i),i);
end
cmf=bsxfun(@minus,csw,min(csw))/nsampl;
cummassflux=bsxfun(@minus,csw,mid)/nsampl;
plot(sw,cummassflux)

cm=b2rcolormap(9);
cm=cm([1:3 6:8],:);
zzi=round(9:9.5:58);
clf;
subplot(2,1,1)
hold on;
for zi=6:-1:1
    plot(sw(:,zzi(zi)),cmf(:,zzi(zi)),'color',cm(zi,:),'linewidth',1.4)
end
legend(sprintf('%3.0f m',h(zzi(6))),sprintf('%3.0f m',h(zzi(5))),sprintf('%3.0f m',h(zzi(4))),...
       sprintf('%3.0f m',h(zzi(3))),sprintf('%3.0f m',h(zzi(2))),sprintf('%3.0f m',h(zzi(1)))    );
set(gca,'fontsize',14,'xlim',[-3 3],'xtick',-3:0.5:3)
title('cumulative volume flux (m s^{-1})')
xlabel('vertical velocity (m s^{-1})')
saveas(gcf,'cu_volflux.eps','epsc')

M=bsxfun(@times,acount,wedges+0.25);

fh=figure;
clf
bar(wedges,M(:,49),'histc');
saveas(fh,'test.png','png')

