function [Vgbin,Ustdbin,nbin]=Pinsky_makelookup_cloudrelative_thickness(V,Z,Zbin,h,htop,hthick,thickbin)
% [Vg,W,Vgbin,Ustdbin,a]=Pinsky_makelookup_cloudrelative(V,Z,Zbin,rho,h,htop,hthick)
% Make the Pinsky et al. (2010) statistical settling velocity lookup table, on height
% coordinate relative to cloud top.
% 
% INPUTS
% V     Doppler velocity
% Z     Reflectivity
% Zbin  reflectivity bins for computing mean and variance of V [default -40:0.25:10]
% 
% OUTPUTS
% Vgbin     Z-dependent average Vg, Pinsky's phi(Zbin)
% Ustdbin   Z-dependent standard deviation of residual velocity U=V-Vgbin(Z)
%           Pinsky's theta(Zbin)
% NOTES
% Bin averages are used rather than 6th order polynomial fit because it is not statistically stable.

% Simon de Szoeke 2013 January 30

if length(Zbin)==1
    Zbin=(-40:Zbin:10)';
end
if length(thickbin)==1
    thickbin=(0:thickbin:1400)';
end

% radar range gate index of cloud top
itop=interp1(h,1:length(h),htop,'nearest');
jtop=91;
shift=-90:10;
nhshift=length(shift);
nh=size(Z,2);

% mean and variance by reflectivity bin
[Vgbin,Ustdbin]=deal(NaN(length(Zbin),length(thickbin),nhshift));
% Ustdbin is Pinsky's theta(Z)
nbin=zeros(length(Zbin),length(thickbin),nhshift);

% cap inputs with NaNs vertically
Z(:,[1 end])=NaN; 
V(:,[1 end])=NaN;

Zshift=NaN(size(Z,1),nhshift);
Vshift=NaN(size(Z,1),nhshift);

% indices for shifted (target) array, 1:101, 91=cloud top
%jj=jtop+shift;
% indices for input (source) array
ii=max(1,min(nh,bsxfun(@plus,itop,shift)));
% bsxfun does addition of column and row vectors with matix mult. rules

% shift Z,V so the cloud top is always at index 91
for ti=1:size(Z,1)
    Zshift(ti,:)=Z(ti,ii(ti,:));
    Vshift(ti,:)=V(ti,ii(ti,:));
end

for hi=1:nhshift
    [tmp1,tmp2,~,tmp3]=binavg2d(Zbin,Zshift(:,hi),thickbin,hthick,Vshift(:,hi),1);
    Vgbin(:,:,hi)=tmp1;
    Ustdbin(:,:,hi)=tmp2;
    nbin(:,:,hi)=tmp3;
end
return