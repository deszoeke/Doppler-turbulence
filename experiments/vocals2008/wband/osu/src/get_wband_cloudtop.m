function [flag,mean_top_bin,median_top_bin,val15_top_bin,val85_top_bin,cloud_fraction]=...
    get_wband_cloudtop(Z,H,noisefloor,time_offset,dt,range_gate);
% [flag,topmean,topmedian,top15,top85,cloud_fraction]=get_wband_cloudtop(Z,H,noisefloor,time_offset,dt);
%
% 2013-03-15 :: VOCALS 2008 :: Simon de Szoeke
%
% Find the W-band radar cloud top in time intervals of dt from
% reflectivity Z>noisefloor located on the vertical coordinate matrix H.
%
% H is the size of Z.
%
% flag values
% 0 no errors detected, cloud top computed.
% 1 cloud fraction=0, cloud top set to NaN.
% 2 cloud fraction>=0 but found no 3 cloud returns with consistent height,
%   possible glitch, provisional cloud top computed.


%cd /Users/sdeszoek/Data/cruises/VOCALS_2008/RHB/Scientific_analysis/programs/VOCALS2008_programs/wband
%addpath('~/Data/cruises/VOCALS_2008/RHB/Scientific_analysis/programs/VOCALS2008_programs/');
read_parameters;

ntimes=size(Z,1);

% Determine cloud top height for each return

% Use a reflectivity criterion input noisefloor that selects clouds
% and eliminates (almost) all noise returns. Further reduce
% nonmeteorological noise and clutter by requiring clouds to be
% vertically homogeneous for 3 range gates.
iscloud1=Z(1:ntimes,:)-repmat(noisefloor,[ntimes 1])>0 & H<2e3 & H>400; % Z criterion for cloud
iscloud2=iscloud1 & circshift(iscloud1,[0 1]); % also require cloud immediately below
iscloud3=iscloud2 & circshift(iscloud1,[0 2]);
iscloud=iscloud3;
isglitch=iscloud1 & ~iscloud3; % 1 and 2-point clouds considered glitches

% cloud top height
Cloudheight=H.*iscloud; % matrix of heights in cloud, 0 out of cloud, still has glitches
soo=sort(Cloudheight(isfinite(Cloudheight)),'descend');
highest=find(abs(soo(1:100)-soo(3:102))<range_gate/2,1,'first'); % highest cloud height for which there are at least 3 retrievals
if ~isempty(highest)
    Cloudheight(Cloudheight>soo(highest))=0; % block too-high glitches
else
    Cloudheight(:)=0; % cloudheight is all glitches
end

% Cloud top height is determined for length dt (60-second) bins:
% The height for which there are at least 3 cloudy returns
% in the minute is the highest allowed top. The cloud top for the
% minute is taken as the mean of highest cloud tops at or below the
% highest allowed top.
bin=0:dt:time_offset(end)+1.0;
nbin=length(bin);
time_ind=floor(interp1(bin,1:nbin,time_offset,'linear'));

% allocate arrays
num_returns=zeros(nbin-1,1);
num_cloud_returns=zeros(nbin-1,1);
num_nonmet_returns=zeros(nbin-1,1);
flag=zeros(nbin-1,1,'int8');
mean_top_bin=NaN+zeros(nbin-1,1);
median_top_bin=mean_top_bin;
std_top_bin=mean_top_bin;
val15_top_bin=mean_top_bin;
val85_top_bin=mean_top_bin;
meanval_top_bin=mean_top_bin;
cldfrac1=mean_top_bin;
cldfrac2=mean_top_bin;
cldfrac3=mean_top_bin;
cloudtop_filealltime=zeros(size(time_offset));

% loop 1 minute bins
for i=1:nbin-1
    in_bin=time_ind==i;
    num_returns(i)=sum(in_bin);
    % is a cloud present?
    num_cloud_returns0=sum(any(iscloud(in_bin,:),2)); % cloud or nonmet
    % Separate cloud-or-nonmeteorological returns into cloud and nonmet by
    % insisting there are 3 returns at given height per minute for clouds.
    if num_cloud_returns0>0
        % cloud height
        C=Cloudheight(in_bin,:);
        soi=sort(C(:),'descend');
        highest=find(abs(soi(1:100)-soi(3:102))<range_gate/2 & soi(1:100)>0,1,'first'); % highest allowed: at least 3 at this level must be found
        if isempty(highest)
            flag(i)=2; % possible clouds, but no cloud with consistent height found
            num_cloud_returns(i)=0;
            num_nonmet_returns(i)=num_cloud_returns0; % assign all to nonmet
        else
            flag(i)=0; % no errors
            num_cloud_returns(i)=num_cloud_returns0-(highest-1);
            num_nonmet_returns(i)=highest-1;
            C(C>soi(highest))=0; % exclude tops above highest allowed
        end
        % proceed to compute cloud top
        cloudtop_filealltime(in_bin)=max(C,[],2);
        valsort=sort(cloudtop_filealltime(in_bin & cloudtop_filealltime>0));
        lensort=sum(in_bin & cloudtop_filealltime>0);
        if lensort>=6 %%needs to be tested
            val15_top_bin(i)=valsort(round(0.15*lensort));
            val85_top_bin(i)=valsort(round(0.85*lensort));
            meanval_top_bin(i)=mean( valsort(round(0.15*lensort):round(0.85*lensort)) );
        else
            val15_top_bin(i)=NaN;
            val85_top_bin(i)=NaN;
            meanval_top_bin(i)=NaN;
        end
        median_top_bin(i)=median(cloudtop_filealltime(in_bin & cloudtop_filealltime>0));
        mean_top_bin(i)=mean(cloudtop_filealltime(in_bin & cloudtop_filealltime>0));
        std_top_bin(i)=std(cloudtop_filealltime(in_bin & cloudtop_filealltime>0));
    else % num_cloud_returns0==0
        flag(i)=1; % cloud fraction=0
        num_cloud_returns(i)=0;
        num_nonmet_returns(i)=0;
        cloudtop_filealltime(in_bin)=NaN;
        median_top_bin(i)=NaN;
        mean_top_bin(i)=NaN;
    end
    if num_cloud_returns<0
        disp('Whoa! Cloudy return number < 0')
        break
    end
    
    % simpler, more liberal cloud fraction counters
    cldfrac1(i)=sum(any(iscloud1(in_bin,:),2))/num_returns(i);
    cldfrac2(i)=sum(any(iscloud2(in_bin,:),2))/num_returns(i);
    cldfrac3(i)=sum(any(iscloud3(in_bin,:),2))/num_returns(i);
    glitchfrac(i)=sum(any(isglitch(in_bin,:),2))/num_returns(i);
    
end % loop 1-min bins
num_clear_returns=num_returns-num_cloud_returns-num_nonmet_returns;
cloud_fraction=num_cloud_returns./(num_cloud_returns+num_clear_returns);




