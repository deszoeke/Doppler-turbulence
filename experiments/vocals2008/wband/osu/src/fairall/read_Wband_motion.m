disp('read wband motion');
vkonsa=[];%ones(10,36000);
vkons=[];
fm = [way_raw_data_wband 'motion_adjT/' year ddds hr '00Kongsberg_adjT.txt'];
FID2 = fopen(fm,'r');
dumx=fgetl(FID2);
dumx=fgetl(FID2);
dumx=fgetl(FID2);
dumx=fgetl(FID2);
for i=1:90%For some reason a fscanf(fid,format,[10,inf]); terminates after 410 lines.  No idea why, but this works. I hate Bill Gates.
    vkonsa=fscanf(FID2,'%2d %*c %2d %*c %4d %2d %*c %2d %*c %g %g %g %g %g',[10,400]);
    vkons=[vkons vkonsa];
    vkonsa=[];
end;
hhk=vkons(4,:)+vkons(5,:)/60+vkons(6,:)/3600;
ii=find(diff(hhk)<=0);
for k=1:length(ii)
    hhk(ii(k)+1)=hhk(ii(k)+1)+1e-6/length(ii)*k;
end;


