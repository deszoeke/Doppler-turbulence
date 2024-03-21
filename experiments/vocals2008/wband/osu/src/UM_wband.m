% DEW KM read in W-Band nc moment data and plot
% range gates are 25 m
% timing is every 0.25 sec
close all
clear
read_parameters_VOCALS

if 0
ddds = input('Input monthday  ','s');
ddd = str2num(ddds);
shrs = input('Input start hr  ','s');
shr = str2num(shrs);
ehrs = input('Input end hr  ','s');
ehr = str2num(ehrs);
end

way_raw_images_wband= 'H:\RHB\radar\wband-UM\Raw_images\';
%ncload('C:\Data\VOCALS\Data\wband\NetCDFRadarmfiles\WbandMom.nc')  % First rain shower

%20082821346MMCRMom.nc  20082821400MMCRMom.nc% 2nd rain shower
way_raw_data_wbandUM = 'H:\RHB\radar\wband-UM\Raw\';
ht = (1:1:120)*.029; % range gates in km
ls([way_raw_data_wbandUM '*1025*.nc'])

%fn = 'H:\RHB\radar\wband-UM\Raw\vocals.94GHz.radar.moments.20081020.150038.218.nc';
fn = 'H:\RHB\radar\wband-UM\Raw\vocals.94GHz.radar.moments.20081025.114500.625.nc'; 
%fn(35:37) = ddds;

for jam = 1:1  %shr:ehr
    if jam<10,
        hr=['0' num2str(jam)];
    else
        hr=num2str(jam);
    end;
 %   fn(38:39) = hr
    ncload(fn)

    tt = time_offset + base_time;
    % calculate time from time_offset and base_time
    for i = 1:length(tt)
        %X = base_time;
        y(i) = tt(i)/(60*60*24);
        z(i) = y(i)+datenum([1970,1,1,0,0,0]);
        dt = datestr(z(i), 31);    %2008-10-07 16:35:53
        year = str2num(dt(1:4));
        month = str2num(dt(6:7));
        day = str2num(dt(9:10));
        hour = str2num(dt(12:13));
        minutes = str2num(dt(15:16));
        secs = str2num(dt(18:19));
        jd(i) = datenum(year,month,day,hour,minutes,secs)-datenum(year-1,12,31);
    end
    jd = (jd-(fix(jd(1))))*24; % subtract off JD and convert to an hour
    Ref = Reflectivity';
    SW = SpectralWidth';
    DV = MeanDopplerVelocity';
    
    if 0
    figure(1)
    %imagesc(jul_day_sonde,[(1:1001)-1]*25/1000,real(theta(1:1001,:)),[280,340]);  % from balloon
    imagesc(jd,ht, Ref(1:120,:),[-60, 20]);  % from balloon
    axis xy
    xlabel('Hrs (UTC)')
    ylabel('Ht (km)')
    title('W-band VOCALS 2008 1635 UTC 7 OCT')
    h = colorbar;
    axes(h)
    ylabel('dBZ')

    figure(2)

    imagesc(MeanDopplerVelocity')
    axis xy
    xlabel('Hrs (UTC)')
    ylabel('Ht (km)')
    title('W-band VOCALS 2008 1635 UTC 7 OCT')
    h = colorbar;
    axes(h)
    ylabel('DopVel (m/s)')

    figure(3)

    imagesc(SpectralWidth')
    axis xy
    xlabel('Hrs (UTC)')
    ylabel('Ht (km)')
    title('W-band VOCALS 2008 1635 UTC 7 OCT')
    h = colorbar;
    axes(h)
    ylabel('SpectralWidth (m/s)')
    end  % skip 3 plots


    figure(4)
    subplot(3,1,1),imagesc(jd,ht, Ref(1:120,:),[-60, 20]);  %
    axis xy
    ylabel('Ht (km)')
    tit = {[cruise num2str(year) ' W-Band'' YD ' ddds ' Hr ' hr]};
    title(tit)
    h = colorbar;
    axes(h)
    ylabel('dBZ')

    subplot(3,1,2),imagesc(jd,ht, DV(1:120,:),[-7, 7]);  %
    axis xy
    ylabel('Ht (km)')
    h = colorbar;
    axes(h)
    ylabel('DopVel (m/s)')


    subplot(3,1,3),imagesc(jd,ht, SW(1:120,:),[0, 5]);  %
    axis xy
    ylabel('Ht (km)')
    xlabel('Hrs (UTC)')
    h = colorbar;
    axes(h)
    ylabel('SpecWid (m/s)')
    print('-djpeg90 ',[way_raw_images_wband '\Wband_' num2str(ddd), '_' hr '.jpg']);

    clear SW DV Ref tt jd
    
end % hr loop
