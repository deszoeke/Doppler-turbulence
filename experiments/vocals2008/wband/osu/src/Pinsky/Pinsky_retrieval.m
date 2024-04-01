function [Vg,W,Vgbin,Ustdbin,Ustd,a]=Pinsky_retrieval(V,Z,Zbin,rho)
% [Vg,W,Vgbin,Ustdbin,a]=Pinsky_retrieval(V,Z,Zbin,rho)
% perform the Pinsky et al. (2010) statistical velocity retrieval from radar Doppler velocity
% 
% INPUTS
% V     Doppler velocity
% Z     Reflectivity
% Zbin  reflectivity bins for computing mean and variance of V [default -40:0.25:10]
% rho   specified height-dependent correlation of W and Vg'    [default 0]
% 
% OUTPUTS
% Vg    Doppler estimate of mean settling velocity
% W     air velocity retrieval
% Vgbin     Z-dependent average Vg, Pinsky's phi(Zbin)
% Ustdbin   Z-dependent standard deviation of residual velocity U=V-Vgbin(Z)
%           Pinsky's theta(Zbin)
% NOTES
% Vgprime=V-Vg-W
% Bin averages are used rather than 6th order polynomial fit because it is not statistically stable.

% Simon de Szoeke 2012 December 28

% default correlation rho and reflectivity bins
if nargin<4
    rho=zeros(1,size(V,2));
    if nargin<3
        Zbin=(-40:0.25:10)';
    end
end  
if length(Zbin)==1
    Zbin=(-40:Zbin:10)';
end
rho=rho(:)'; % row vector

nh=size(Z,2);

% mean and variance by reflectivity bin
[Vgbin,Ustdbin]=deal(NaN(length(Zbin),nh));
% Ustdbin is Pinsky's theta(Z)

nbin=zeros(length(Zbin),nh);
indx=zeros(size(Z));
[Vg,Ustd]=deal(NaN(size(Z)));
for hi=2:100
    [Vgbin(:,hi),Ustdbin(:,hi),~,nbin(:,hi)]=binavg(Zbin,Z(:,hi),V(:,hi),8);
    
%     % 6th order polynomial coefficient
%     isf=isfinite(Vgbin(:,hi));
%     Cphi(:,hi)=polyfit(Zbin(isf)'+0.125,Vgbin(isf,hi),6);
%     Ctheta(:,hi)=polyfit(Zbin(isf)'+0.125,Ustdbin(isf,hi),6);
    
    % skip the polyfits and use Vgbin and Ustdbin directly
    % lookup Vg and Ustd=theta(Z(h,t)) for all h,t
    indx=floor(interp1(Zbin,1:length(Zbin),Z(:,hi),'linear'));
    indx(isnan(indx))=length(Zbin);
    indx(indx==0)=length(Zbin);
    Vg(:,hi)=min(0,Vgbin(indx,hi)); % enforce Vg<=0
    Ustd(:,hi)=Ustdbin(indx,hi);
end
Vg(isnan(Ustd))=NaN;
% consider option to attribute some of this Vg to Wbar(Z).

% residual velocity
U=V-Vg;

% linear weighted attribution of residual to W and Vg'
S1=nanmean(1./Ustd);
S2=nanmean(1./(Ustd.*Ustd));
% for correlation rho=<WVg'>=0
% a0=S1./S2;

% general case for zero or nonzero correlation rho
a0=S1./S2-rho./S2.*sqrt((S2-S1.*S1)./(1-rho.*rho));
a=repmat(a0,[size(V,1),1])./Ustd;
% separate W and Vgprime from U
W=a.*U;
Vgprime=(1-a).*U;

% Vgprime=U-W=V-Vg-W;
return