function [Vg,W,Ustd,a]=Pinsky_retrieval_cloudrelative(V,Z,Zbin,Vgbin,Ustdbin,nbin,rho,h,htop,foot)
% [Vg,W,Vgbin,Ustdbin,a]=Pinsky_retrieval(V,Z,Zbin,Ustdbin,rho,h,htop,hthick)
% perform the Pinsky et al. (2010) statistical velocity retrieval from radar Doppler velocity
% 
% INPUTS
% V     Doppler velocity
% Z     Reflectivity
% Zbin  reflectivity bins for computing mean and variance of V [default -40:0.25:10]
% Vgbin   lookup table for Vg(Z,h-htop)=phi(Z)
% Ustdbin   lookup table for theta(Z,h-htop)=<U'^2>_Z^1/2
% nbin  number of observations in heat h,Z bin
% rho   specified height-dependent correlation of W and Vg'    [default 0]
% h     range height (m)
% htop  cloud top height (m)
% foot  minimum Ustd (m/s) averaged to get S1,S2, avoids numerical glitches
%
% OUTPUTS
% Vg    Doppler estimate of mean settling velocity
% W     air velocity retrieval
% Ustd  standard deviation of residual velocity U=V-Vgbin(Z)
%
% NOTES
% Vgprime=V-Vg-W
% Bin averages are used rather than 6th order polynomial fit because the fit is not statistically stable.

% Simon de Szoeke 2012 December 28

rho=rho(:)'; % row vector

% compute bin average a0 from the height-relative table
sigma=Ustdbin;
sigma(sigma==0 | nbin<12)=NaN; % avoid divide overflow, preserve structure of sigma
s1=nansum(nansum(nbin./sigma,3),1)./nansum(sum(nbin,3),1);
s2=nansum(nansum(nbin./(sigma.*sigma),3),1)./sum(sum(nbin,3),1);
a0=s1./s2; % function of top-relative-height (already in shifted coordinate)

% radar range gate index of cloud top
nh=length(h);
iihtop=isfinite(htop);
itop=ones(size(htop));
itop(iihtop)=interp1(h,1:nh,htop(iihtop),'nearest');
itop(~iihtop)=1; % kluge for undefined cloud tops
%jtop=91; % the index of the cloud top once shifted
shift=-90:10;
nhshift=length(shift);

% mean and variance by reflectivity bin
% Ustdbin is Pinsky's theta(Z)
% nbin=zeros(length(Zbin),nh);
% indx=zeros(size(Z));
Zshift=zeros(size(Z,1),nhshift);

% indices for input (source) array
ii=max(1,min(nh,bsxfun(@plus,itop,shift)));
% bsxfun does addition of column and row vectors with matix mult. rules

% shift Z,V so the cloud top is always at index 91
for ti=1:size(Z,1)
    Zshift(ti,:)=Z(ti,ii(ti,:));
%     Vshift(ti,:)=V(ti,ii(ti,:));
end

[Vgshift,Ustdshift,Wbarupshift]=deal(NaN(size(Zshift)));
[Vg,Ustd,Wbarup,a]=deal(NaN(size(Z)));
for hi=1:nhshift
    % lookup Vg and Ustd=theta(Z(h,t)) for all hrel,t [cloud-top rel. coordinate]
    indx=floor(interp1(Zbin,1:length(Zbin),Zshift(:,hi),'linear'));
    indx(isnan(indx))=length(Zbin);
    indx(indx==0)=length(Zbin);
    Vgshift(:,hi)=min(0,Vgbin(indx,hi));     % enforce Vg<=0
    Wbarupshift(:,hi)=max(0,Vgbin(indx,hi)); % mean upward "fall" velocity
    Ustdshift(:,hi)=Ustdbin(indx,hi);
end
Vgshift(isnan(Ustdshift))=NaN;
% consider attributing some of this Vg to Wbar(Z).

% % linear weighted attribution of residual to W and Vg'
% % compute a0shift(hrel) in the cloud relative coordinate
% Ustdtmp=max(foot,Ustdshift);
% Ustdtmp(isnan(Ustdshift))=NaN;
% S1shift=nanmean(1./Ustdtmp);
% S2shift=nanmean(1./(Ustdtmp.*Ustdtmp));
% % general case for zero or nonzero correlation rho
% %a0=S1./S2-rho./S2.*sqrt((S2-S1.*S1)./(1-rho.*rho));
% a0shift=S1shift./S2shift;

% now a0 comes from lookup table above
ashift=repmat(a0,[size(V,1),1])./max(foot,Ustdshift);

% shift back to standard height coordinate using ii
for ti=1:size(Z,1)
    Vg(ti,ii(ti,:))=Vgshift(ti,:);
    Ustd(ti,ii(ti,:))=Ustdshift(ti,:);
    Wbarup(ti,ii(ti,:))=Wbarupshift(ti,:);
    a(ti,ii(ti,:))=ashift(ti,:);
end

% residual velocity
U=V-Vg;

% doing second stage separation in ground-relative h works, but unintended
% SPdeS moved up to cloud-relative 2013-04-04

% separate W and Vgprime from U
W=a.*U; % SPdeS Apr4 limit a to [0 1]
%W=min(1,a).*U; % SPdeS Apr4 limit a to [0 1]
% Vgprime=(1-a).*U;

% Vgprime=U-W=V-Vg-W;
return