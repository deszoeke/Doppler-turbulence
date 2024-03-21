% plot_dBZ.m
% 
% Plot theoretical dBZ contours for drop concentration and radius.
% This will depend on the drop size distribution--try:
% 1. Dirac delta radius distribution at <D>
% 2. Lognormal radius distribution
%
% Simon de Szoeke

% 1. Dirac delta diameter distribution at <D>
D=[1:.5:9 10:5:90 100:50:900 1e3]'; % micron
N=10.^(1:.025:3.5); % cm-3
Z=(D*1e-3).^6 * N*1e6; % D:micron-->mm, N:cm-3-->m-3
dBZ1=10*log10(Z);

% [c,h]=contour(D,N,dBZ1');
% set(gca,'xscale','log','yscale','log')
% xlabel('diameter (\mum)')
% ylabel('number (cm^{-3})')
% title({'delta distribution' 'dBZ'})
% clabel(c,h,'manual')

% 2. Lognormal radius distribution
% assume D is the mode diameter D_0
sigx=0.35;
Z=(D*1e-3).^6 * N*1e6 * exp(18*sigx^2);
dBZ2=10*log10(Z);

[c,h]=contour(D,N,dBZ2');
set(gca,'xscale','log','yscale','log')
xlabel('diameter (\mum)')
ylabel('number (cm^{-3})')
title({'lognormal distribution, \sigma_x=0.35' 'dBZ'})
clabel(c,h,'manual')
box on

print -deps plot_dBZ_lnorm.eps
print -dpng plot_dBZ_lnorm.png

% lognormal distribution of width 0.35 accounts for 10 dBZ higher reflectivity than delta distribution (width=0).
x=0:0.05:0.75;
DdBZ = 78.17*x.^2;

plot(x,DdBZ,'linewidth',1.4)
axis tight; grid on
hold on
plot(.35,78.17*.35^2,'x','markersize',10)
ylabel('\Delta dBZ')
xlabel('logarithmic width \sigma_x')
print -depsc dBZ_logwidth.eps