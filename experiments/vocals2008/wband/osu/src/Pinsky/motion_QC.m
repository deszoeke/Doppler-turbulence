motion_status
% 0 OK
% 1 no data
% 2 bias (<1 degree)
% 3 noisy
% 4 motion adj. off or failed

% Doppler Velocity shifts AND Suspect sea motion contamination
x1=[[319 18
320 05
320 11
320 14
322 02
322 05
322 14
329 13.88
329 14.25
330 10]];

% Suspect ~10s sea motion contamination
x2=[[324 15
324 18
326 08
326 11
327 23
330 03
331 09]];

% adjusted Doppler velocity is missing, usually around 15-18 minutes of the hour
x3=[[328 02
328 05
328 08
328 11
328 14
330 05
330 08
330 11
330 14
335 09]];

clf
imagesc(statushour,statusday,Status)
caxis([-.5 4.5])
colormap([1 1 1; 0 0 0; 0 .3 .8; .9 0 0; .7 .6 0])
cb=colorbar;
set(cb,'ytick',0:4,'yticklabel',{'OK','no data','bias','noisy','failed/off'})
set(gca,'fontsize',14,'xtick',0:4:24)
xlabel('UTC hour')
ylabel('2008 yearday')
title('VOCALS W-band motion status & Doppler velocity notes')
hold on
plot(x1(1:end-2,2),x1(1:end-2,1),'gv','markersize',10)
plot(x1(end-1:end,2),x1(end-1:end,1),'g^','markersize',10)
plot(x2(:,2),x2(:,1),'gs','markersize',10)
plot(x3(:,2),x3(:,1),'kx','markersize',10)
legend('V shift down','V shift up','suspect 10s sea contamination','V missing',0)

print -dpng motion_QC.png