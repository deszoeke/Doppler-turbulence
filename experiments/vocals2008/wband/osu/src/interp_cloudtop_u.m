% NOT SHIP RELATIVE! 

%stday=281;
%enday=336;

% return Ui, yday_midhour
ydaysond=ncread([way_proc_data_balloon 'VOCALS2008_soundings_z_v4.2.nc'],'yday');
zsond=ncread([way_proc_data_balloon 'VOCALS2008_soundings_z_v4.2.nc'],'height');
% interp scalar wind speed for dissipation calc.
wndspd=ncread([way_proc_data_balloon 'VOCALS2008_soundings_z_v4.2.nc'],'wndspd');
yday_midhour=(stday:1/24:enday)'+1/48;
xi=interp2(ydaysond(:),zsond(:),wndspd,yday_midhour',h);
Ui=xi(:,1:end-1);

%That was easy!

% u=ncread([way_proc_data_balloon 'VOCALS2008_soundings_z_v4.2.nc'],'u');
% v=ncread([way_proc_data_balloon 'VOCALS2008_soundings_z_v4.2.nc'],'v');
% % diurnal cycles in u,v
% subplot(2,1,1)
% pcolor(ydaysond,zsond(1:300),u(1:300,:)); shading flat; hold on;
% plot(ydaysond,1500+200*u(200,:))
% subplot(2,1,2)
% pcolor(ydaysond,zsond(1:300),v(1:300,:)); shading flat; hold on;
% plot(ydaysond,1500+500*v(200,:))