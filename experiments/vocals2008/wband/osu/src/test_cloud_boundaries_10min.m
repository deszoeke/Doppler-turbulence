ncpath='/Users/sdeszoek/Data/cruises/VOCALS_2008/RHB/radar/wband/Processed/cloudheight/';
% note netcdflib does not accept ~/ as a path
ncf=[ncpath 'VOCALS2008CloudBoundaries10min_0_2.nc'];

yday=nc_varget(ncf,'yday');
lat=nc_varget(ncf,'lat');
lon=nc_varget(ncf,'lon');
cloudtop=nc_varget(ncf,'cloudtop');
sondecloudtop=nc_varget(ncf,'sondecloudtop');
sondeinterpflag=nc_varget(ncf,'sondeinterpflag');
cloudbase=nc_varget(ncf,'cloudbase');
ncloud=nc_varget(ncf,'numceilcloud');

plot(lon,lat); hold on
isonde=find(sondeinterpflag==1);
for iso=1:length(isonde) 
    hm(iso)=rectangle('Curvature',[1 1],'Position',...
        [lon(isonde(iso))-0.5*sondecloudtop(isonde(iso))/1e3 lat(isonde(iso))-0.5*sondecloudtop(isonde(iso))/1e3...
         [0 0]+sondecloudtop(isonde(iso))/1e3]);
end
axis equal

clf
plot(yday,sondecloudtop,'r.')
hold on
plot(yday,sondecloudtop,'r-')
plot(yday,cloudbase,'b.')
plot(yday,cloudbase,'b-')
plot(yday,cloudtop,'.','color',0.7*[1 1 1])
plot(yday,cloudtop,'-','color',0.7*[1 1 1])


%rain contamination is a problem for cloud base retrieval