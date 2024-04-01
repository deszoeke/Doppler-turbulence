function [time,w,pitch,roll]=read_kongsberg(filename)
% [time,w]=read_kongsberg(filename)
% Reads the file filename and returns the decimal yearday time
% and vertical velocity w.
%
% (c) 2009-04-22 Simon de Szoeke

%filename=[way_raw_data_wband 'motion/20083301247Kongsberg.txt']
%filename='/Users/sdeszoek/Data/cruises/VOCALS_2008/RHB/radar/wband/Raw/motion_adjT/20083300600Kongsberg_adjT.txt'
fid=fopen(filename);
%frewind(fid);
adjT=strcmp(filename(end-7:end),'adjT.txt');

% read header
timestamp=fscanf(fid,'Kongsberg output log starting at %4s%3s%2s%2s  (YYYYJJJHHMM)\n',4);
if adjT
    fscanf(fid,'%*s',11); % longer header
else
    fscanf(fid,'%*s',9);
end

if adjT
    A=fscanf(fid,'%*2d/%*2d/%*4d  %2d:%2d:%6f %11e %11e %11e %11e %11e',[8,inf])';
else
    A=fscanf(fid,'%*2d/%*2d/%*4d  %2d:%2d:%6f %10e %10e %10e %10e',[7,inf])';
end
yday=sscanf(timestamp,'%*4c%3d%*4c',3);
hour=A(:,1);
minute=A(:,2);
second=A(:,3);
pitch=A(:,4);
roll=A(:,5);
w=A(:,6);
time=yday+datenum(0,0,0,hour,minute,second);
%time_offset=A(:,8); % computed by Alan Brewer and appended to Kongsberg_adjT.txt files

fclose(fid);
% 10/06/2008  12:47:45.422  1.957e-02      1.618e-02      -2.220e-04      -2.220e-04