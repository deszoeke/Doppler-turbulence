% plot diurnal decoupling figure with time-height of
% log10(dissipation)
% vertical velocity variance
% vertical velocity skewness

cd('/home/farfalla/data1/sdeszoek/VOCALS/RHB/Scientific_analysis/programs/VOCALS2008_programs/wband/Pinsky')

load('radarturb8.mat');
%tt = load('cloudtopturb6.mat');
p1 = load('/home/farfalla/data1/sdeszoek/VOCALS/RHB/Scientific_analysis/programs/VOCALS2008_programs/wband/Pinsky/retrieval/Pinsky_wRetrieval2008_335_21.mat'); % for height coordinate h
h = p1.h;

% mask variance and skewness
mask = zeros(size(countk));
mask(countk<800) = NaN;

% local hour index for diurnal compositing
indx=binall2(0:24,mod(ydayw10*24-5,24)); % local hour index 1:24; 1 corresp to 0-1 local

% synthesis data
sf = '/home/farfalla/data1/sdeszoek/VOCALS/synthesis/stratus2008_10min_micro.nc';
yday = ncread(sf, 'yday');
% read synthesis variables and interpolate to ydayw10
sread = @(s) interp1(yday, ncread(sf, s), ydayw10);
top  = sread('top');
zlcl = sread('zlcl');
zbmd = sread('zbm');

zlcld = dmn( sread('zlcl'), indx );
topd  = dmn( sread('top' ), indx );
zbmd  = dmn( sread('zbm' ), indx );

% align top
sz = size( epsilonk );
% one attempt
ktop = get_ktop( epsilonk );
% synthesis top, which agrees mostly
ksyntop = 0*ydayw10;
isf = isfinite(top);
fisf = find(isf);
for it = 1:length(fisf)
    ksyntop(fisf(it)) = find( h < top(fisf(it)), 1, 'last');
end
%epsaligntop = vshift( epsilonk, sz(2)-ktop );
epsaligntop = vshift( epsilonk, -ktop );
pcolor(ydayw10, h(1:100)-h(98), circshift(log10(epsaligntop),[0, -2])'); shading flat; caxis([-5 -2.5]); set(gca,'color',0.7+[0 0 0])

epsd = dmn(vshift(epsilonk, -ktop), indx);
vard = dmn(vshift(mask+wvar,     -ktop), indx);
skwd = dmn(vshift(mask+wskw,     -ktop), indx);
ktopd = round(dmn( ktop(ktop>0), indx(ktop>0) ));

% test the cloud top height diagnoses
plot(ydayw10(ktop>0), h(ktop(ktop>0)),'.')
hold on
plot(318:(1/24):319, h(ktopd([1:24 1])),'r-') 
plot(ydayw10, sread('top'))


idi = [1:24 1:12];

ntitle = @(a, s) title(a, s, 'fontweight','norm');
b2rcolormap(13);
clf
ax(1)=subplot(3,1,1);
x = vard;
pcolor(0:35, (h(5:104)-h(100))/1e3, circshift(x(idi,:),[0 -4])'); shading flat
%pcolor(0:35, (h(1:100)-h(100))/1e3, circshift(x(idi,:),[0 -4])'); shading flat
colorbar
hold on
plot(0.5:35.5, -(topd(idi)-zbmd(idi))/1e3, 'k')
plot(0.5:35.5, -(topd(idi)/1e3-zlcld(idi)), 'k')
ax(2)=subplot(3,1,2);
x = skwd;
pcolor(0:35, (h(5:104)-h(100))/1e3, circshift(x(idi,:),[0 -4])'); shading flat
caxis([-0.75 0.75])
colorbar
hold on
plot(0.5:35.5, -(topd(idi)-zbmd(idi))/1e3, 'k')
plot(0.5:35.5, -(topd(idi)/1e3-zlcld(idi)), 'k')
ax(3)=subplot(3,1,3);
x = log10(epsd);
pcolor(0:35, (h(5:104)-h(100))/1e3, circshift(x(idi,:),[0 -4])'); shading flat
colorbar
hold on
plot(0.5:35.5, -(topd(idi)-zbmd(idi))/1e3, 'k')
plot(0.5:35.5, -(topd(idi)/1e3-zlcld(idi)), 'k')
set(ax(:), 'tickdir','out', 'color',0.7+[0 0 0], 'ylim',[-2 0.1]) 
ntitle(ax(1), 'vertical velocity variance (m^2 s^{-2})')
ntitle(ax(2), 'vertical velocity skewness')
ntitle(ax(3), 'log_{10} dissipation (m^2 s^{-3})')
ylabel(ax(2), 'height below cloud top (km)')
set(ax(:),'fontsize',14, 'fontname','Helvetica')
xlabel(ax(3), 'local hour')
set(gcf,'inverthardcopy','off', 'color',[1 1 1])
print -depsc diel_decoupling_var_skw_eps.eps

