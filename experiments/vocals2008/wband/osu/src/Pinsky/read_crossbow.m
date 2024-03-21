function [time,w,pitch,roll]=read_crossbow(filename)
% [time,w,pitch,roll]=read_crossbow(filename)
% Reads the file filename and returns the decimal yearday time
% and vertical velocity w.
%
% (c) 2009-04-22 Simon de Szoeke

%filename='/Users/sdeszoek/Data/cruises/VOCALS_2008/RHB/radar/wband/Raw/motion/20083212100Crossbow.txt'
fid=fopen(filename);
%frewind(fid);

% read header
timestamp=fscanf(fid,'Crossbow output log starting at %4s%3s%2s%2s  (YYYYJJJHHMM)\n',4);
fscanf(fid,'%*s',9); % discard column headers

A=fscanf(fid,'%*2d/%*2d/%*4d  %2d:%2d:%6f %e %e %e %e %e %e %e',[10,inf])';
fclose(fid);

yday=sscanf(timestamp,'%*4c%3d%*4c',3);
hour=A(:,1);
minute=A(:,2);
second=A(:,3);
roll=A(:,4);
pitch=A(:,5);
rollrate=A(:,6);
pitchrate=A(:,7);
xacc=A(:,8);
yacc=A(:,9);
zacc=A(:,10);

time=yday+datenum(0,0,0,hour,minute,second);

% compute vertical vel from acceleration
% 1. integrate accelerations [units?]
% 2. mean velocity is 0 -- so highpass filter
margin=3360;
coef=-4.354e-5; % empirical coefficient for getting units in m/s
fs=20; % Hz
% nyquist=fs/2
wf=2/fs/60; % 60s is the period of the cutoff frequency
w=coef*highpass(cumsum(zacc),5,wf);

% don't need bc interpolating to kongsberg, then radar time anyway
% % interpolate over gaps shorter than 3 s
% ii=isfinite(w+time);
% ibad=diff(time)>3/60/60/24;
% wsave=w(ibad);
% w(ibad)=NaN;
% w(~ii)=interp1(time(ii),w(ii),time(~ii));
% w(ibad)=wsave;

w([1:margin end-margin+1:end])=NaN;

return

if 0
%% compare to kongsberg
filename='/Users/sdeszoek/Data/cruises/VOCALS_2008/RHB/radar/wband/Raw/motion_adjT/20083212100Kongsberg_adjT.txt'
[ktime,kw]=read_kongsberg(filename);

plot(ktime,kw)
hold on
plot(time,-4.5e-5*w)

cw=interp1(time,w,ktime);
ii=isfinite(kw+cw);
P=polyfit(kw(ii),cw(ii),2);
wship=-4.354e-5*w;
end