% This is for zprof files written in the field
%--------------------------------------------------------------------------
%	readZProfProfilesTxt.m
%	This function reads in the zprof files (from a chosen directory)
%--------------------------------------------------------------------------
%     z_alt				- altitude (ASL)
%     atmosVarNsSpec	- atmospheric variance estimated using spectrum
%     atmosVarXCorr		- atmospheric variance estimated using xcorr/ACF
%     ACFstdv			- standard deviation estimated using xcorr/ACF
%     nsSpecStdv		- standard deviation estimated using spectrum
%     totalVar			- total variance as estimated using the ACF lag0
%	  mnSNR				- mean SNR for each range gate (within thresholds)
%     hour				- average hour of the zenith file
%     minute			- and so on
%     sec
%     month
%     day
%     year
%     alt
%     lat
%     lon
%     JulDay

function [out_zProf] = readZProfProfilesTxt(inDir,startDateNum,stopDateNum)
% inDir = 'H:\CSD3\2008_VOCALS\field results\textprofiles\';
% startDateNum = datenum([2008,11,27,0,0,0]);	% change as appropriate [yyyy,mm,dd,hh,MM,ss]
% stopDateNum = datenum([2008,11,27,12,0,0]);	% change as appropriate
fileTypeStr = 'zprof';
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
			numrows(jk)		= X(1); 
			hour(jk)		= X(2);
			minute(jk)		= X(3);
			sec(jk)			= X(4);
			month(jk)		= X(5);
			day(jk)			= X(6);
			year(jk)		= X(7);
			alt(jk)			= X(8);
			lat(jk)			= X(9);
			lon(jk)			= X(10);
			junk			= fgetl(fid);  % read header
			clear datrow;
			if numrows(jk) == 0; continue; end;
			for km = 1:numrows(jk)
				datrow(km,:) = str2num(fgetl(fid)); % getting the data rows for this section
			end			
			indSet = indx:indx+numrows(jk)-1;
			indSet_ab(jk,:) = [indx,indx+numrows(jk)-1];
            
            z_alt   		= datrow(:,1)';
			atmosVarNsSpec  = datrow(:,2)';
			atmosVarXCorr 	= datrow(:,3)';
			ACFstdv 		= datrow(:,4)';
			nsSpecStdv  	= datrow(:,5)';
			totalVar 		= datrow(:,6)';
			mnSNR  		 	= datrow(:,7)';
		
            
        			%load into profdata structure
			num_actual_profs = num_actual_profs +1;
			nap = num_actual_profs;
            zProf(nap).decTime = datenum([year(jk),month(jk),day(jk),hour(jk),minute(jk),sec(jk)]);
            zProf(nap).altitude = alt(jk);
            zProf(nap).latitude = lat(jk);
            zProf(nap).longitude = lon(jk);
            zProf(nap).height = z_alt;
            zProf(nap).atmosVarNsSpec = atmosVarNsSpec;
            zProf(nap).atmosVarXCorr = atmosVarXCorr;
            zProf(nap).ACFstdv = ACFstdv;
            zProf(nap).nsSpecStdv = nsSpecStdv;
            zProf(nap).totalVar = totalVar;
            zProf(nap).mnSNR = mnSNR;
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

%now cut it down to the correct profiles within what we have read
nap = 0;
for jk = 1 : numRecs
	if zProf(jk).decTime < stopDateNum & zProf(jk).decTime > startDateNum 
		nap = nap + 1;
		out_zProf(nap) = zProf(jk);
	end
end
if ~exist('out_zProf') ;
    out_zProf = [];
end;

