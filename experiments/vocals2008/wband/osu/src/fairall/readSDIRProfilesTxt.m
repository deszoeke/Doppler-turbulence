%--------------------------------------------------------------------------
%	readSDIRProfilesTxt.m
%	This function reads in the spdir files (from a chosen directory)
%--------------------------------------------------------------------------

function [out_profData] = readSDIRProfilesTxt(inDir,startDateNum,stopDateNum)
% clear
% sondeDir = 'C:\Data\';
% inDir = 'D:\Data\HRDL\text_profiles\windprofs\'
% startDateNum = datenum([2006,8,4,0,0,0]);	% change as appropriate [yyyy,mm,dd,hh,MM,ss]
% stopDateNum = datenum([2006,8,4,12,0,0]);	% change as appropriate
fileTypeStr = 'spdir';
fprintf('Time Period : %s through %s \n',datestr(startDateNum,'HH:MM mm/dd/yy'),datestr(stopDateNum,'HH:MM mm/dd/yy'));

xdir = cd;
fns = dir(fullfile(inDir,'*.txt'));
nfiles = length(fns);
num_actual_profs = 0;
num_act_files = 0;
for jk = 1:nfiles
	[fl,rem] = strtok(fns(jk).name,'_');
	[fileDateStr,rem] = strtok(rem,'_');
	[fl,rem] = strtok(rem,'.');
	if ~strcmp(strtok(fl,'_'),fileTypeStr)
		continue;		%this is not a legit file tpye
	end
	[fl,rem] = strtok(rem,'.');
	[startTimeStr,rem] = strtok(fl,'to');
	[stopTimeStr,rem] = strtok(rem,'to');

	num_act_files = num_act_files + 1;
	filenms{num_act_files} = fns(jk).name;
	startFileDateNum(num_act_files) = dateNum(strcat(fileDateStr,startTimeStr),'yyyymmddHHMM');
	stopFileDateNum(num_act_files) = dateNum(strcat(fileDateStr,stopTimeStr),'yyyymmddHHMM');
    
%     fprintf('%s   %s - %s\n',fns(jk).name,datestr(startFileDateNum(num_act_files)),datestr(stopFileDateNum(num_act_files)));
end
if ~exist('filenms')
	fprintf('No files of the type %s in %s\n',fileTypeStr,inDir);
	out_profData = [];
	return;
end
[Y,yi] = sort(startFileDateNum);
filenames = filenms(yi)';
startFDN = startFileDateNum(yi); % this is to pick out the files
stopFDN = stopFileDateNum(yi);
gddex = find(startFDN < stopDateNum & stopFDN > startDateNum); 
nfiles = length(gddex );
if  nfiles == 0
	fprintf('No files meet the time requirement\n');
	out_profData = [];
	return
end
filenames = filenames(gddex);
jk = 1;
indx = 1;

for mk = 1:nfiles
	fname = fullfile(inDir,filenames{mk});
	fid = fopen(fname);
	if fid == -1
		fprintf('\Unable to open: %s. \n\n',fname)
		out_profData = [];
		return
	else
		fprintf('Reading data from file %s.\n',filenames{mk});
	end
	myeof = 0;
	while ~feof(fid) & myeof == 0
		str1 = fgetl(fid);				% at start of each section
		if str1 == -1
			myeof = 1;
		else
			%numrows hr min sec mon dy yr el std_thr filenum
			X = str2num(str1);
			if X(1) == -1		%empty profile
				str1 = fgetl(fid);	%read one more line to line up on the next profile
				continue;
			end
			numrows(jk)		= X(1)+1; 
			hour(jk)		= X(2);
			minute(jk)		= X(3);
			sec(jk)			= X(4);
			month(jk)		= X(5);
			day(jk)			= X(6);
			year(jk)		= X(7);
			el(jk)			= X(8);
			alt(jk)			= X(9);
			lat(jk)			= X(10);
			lon(jk)			= X(11);
			a(jk)			= X(12);
			b(jk)			= X(13);
			junk			= fgetl(fid); % read header line
			clear datrow;
			
			for km = 1:numrows(jk)
				datrow(km,:) = str2num(fgetl(fid)); % getting the data rows for this section
			end			
			indSet = indx:indx+numrows(jk)-1;
			indSet_ab(jk,:) = [indx,indx+numrows(jk)-1];
            
            z_alt   		 	= datrow(:,1)';
			windSpeed		 	= datrow(:,2)';
			windDir			 	= datrow(:,3)';
			windStd			 	= datrow(:,4)';
            [r,c] =size(datrow);
            if c == 6  % we have added two additional rows - stdAmp and stdAng
                ampStd			 	= datrow(:,5)';
                angStd			 	= datrow(:,6)';
            else
                ampStd			 	= nan(1,length(windDir));
                angStd			 	= nan(1,length(windDir));
            end

			%load into profdata structure
			num_actual_profs = num_actual_profs +1;
			nap = num_actual_profs;
			profData(nap).z_alt = z_alt;
            profData(nap).dAmp = ampStd;
            profData(nap).dAngle = angStd;
            profData(nap).windSpeed = windSpeed;
            profData(nap).windDir = windDir;
            profData(nap).windStd = windStd;
            profData(nap).numrows = numrows(jk);
            profData(nap).hour = hour(jk);
			profData(nap).minute = minute(jk);
			profData(nap).sec =	sec(jk);
			profData(nap).month = month(jk);
			profData(nap).day = day(jk);
			profData(nap).year = year(jk);
			profData(nap).beamEl = el(jk);
			profData(nap).alt = alt(jk);
			profData(nap).lat = lat(jk);
			profData(nap).lon = lon(jk);
			profData(nap).a = a(jk);
			profData(nap).b = b(jk);
			profData(nap).JulDay = datenum([profData(nap).year,profData(nap).month, ...
									   profData(nap).day, profData(nap).hour,...
									   profData(nap).minute,profData(nap).sec]);
			if length(find(~isnan(z_alt)))==0
				fprintf('no values for this file. jk: %i\n',jk);
			else
				jk = jk+1;
			end
		end
	end
	fclose(fid);

end
numRecs = jk-1;
axmax = 4000;

%now cut it down to the correct profiles within what we have read
nap = 0;
for jk = 1 : numRecs
	if profData(jk).JulDay < stopDateNum & profData(jk).JulDay > startDateNum 
		nap = nap + 1;
		out_profData(nap) = profData(jk);
	end
end
if ~exist('out_profData') ;
    out_profData = [];
end;
return

%--------------------------------------------------------------------------
% now load the sonde data


sondfname = fullfile(sondeDir,sondFile);
fid = fopen(sondfname);
if fid == -1
	fprintf('\Unable to open: %s. \n\n',sondfname)
    return
else
	fprintf('Reading sonde data from file %s.\n',sondFile);
end

for mk = 1:16,	junk = fgetl(fid); end

jk = 1;
myeof = 0;
while ~feof(fid) & myeof == 0
	str1 = fgetl(fid);					% at start of each section
	if str1 == -1
		myeof = 1;
	else
		X = str2num(str1);
		time_sec(jk) = X(1); 
		press_mb(jk) = X(2);
		temp_C(jk) = X(3);
		RH_pct(jk) = X(4);
		Wdir_deg(jk) = X(5);
		Wspd_mps(jk) = X(6);
		dz_mps(jk) = X(7);
		Lon_deg(jk) = X(8);
		Lat_deg(jk) = X(9);
		Alt_m(jk) = X(10);
		Sa_nums = fgetl(fid);
		jk = jk+1;
	end
end
fclose(fid);

numSondeRecs = jk-1;
axmax = 4000;
badInds1 = find(Wspd_mps == -999 | Alt_m == -999);
badInds2 = find(Wdir_deg == -999 | Alt_m == -999);
Wspd_mps(badInds1) = nan;
Wdir_deg(badInds2) = nan;
Alt_m(badInds1) = nan;

%--------------------------------------------------------------------------
% plot both data sets together
%
figure(1); clf;
for jk = 1:numRecs
%	jk
	hrstr = num2str(hour(jk),2);
	minstr = num2str(minute(jk),2);
	secstr = num2str(sec(jk),2);
	if length(minstr) == 1; minstr = strcat('0',minstr);end
	if length(secstr) == 1; secstr = strcat('0',secstr);end
	titlestr = strcat(hrstr,':',minstr,':',secstr);
	stinds = indSet_ab(jk,1):indSet_ab(jk,2);
	subplot(121);
	plot(Wspd_mps,Alt_m,windSpeed(stinds),z_alt(stinds),'bd');
	hold on;
	for lk = stinds(1):stinds(end)
		xdims = windSpeed(lk)+[-windVar(lk)/2,windVar(lk)/2];
		ydims = [z_alt(lk),z_alt(lk)];
		hp=plot(xdims,ydims,'r');
		set(hp,'linewidth',2);
	end
	hold off
	axis([0,27,0,axmax]);
	title(titlestr);
	xlabel('Wind Speed (m/s)');
	ylabel('Altitude (m)');
	subplot(122);
	plot(Wdir_deg,Alt_m,windDir(stinds),z_alt(stinds),'bd'); 
	axis([0,360,0,axmax]);
	title(sondFile)
	xlabel('Wind Direction (deg.)');
	pause
end

