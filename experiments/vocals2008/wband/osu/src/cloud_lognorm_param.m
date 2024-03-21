function [SigX, ModalRadius, N]=cloud_lognorm_param(DopplerV, DopplerVWidth, Z)
% [SigX, ModalRadius, N]=cloud_lognorm_param(DopplerV, DopplerVWidth, Z);
%  ln(m)   m                                        m/s            m/s          m^3
%
% Find the parameters SigX. ModalRadius, and N for lognormally distributed
% drizzle drops,
%
% n_cloud(x) = N/(SigX*sqrt(2*pi))*exp[-(x-x0)^2/2*SigX^2]
%
% with x=ln(r), x0=ln(ModalRadius), from a single realization Doppler vertical drop
% velocity and Doppler velocity width. The Stokes fall-speed radius relation for
% cloud drops is adopted Vfall=c*r, c=2*g*rhow/(9*visc_dyn).
% The Doppler radar velocity is actually the sum of the air w and the fall
% speed. The air velocity and air velocity width are assumed to be zero in this
% function, so Doppler velocities are assumed to indicate only fall velocity.
% Z and N are optional.
% 
% N is unreliable from this method, but could be compared  to LWC from microwave.
% This retrieval unlikely to be right without quantifying and including contribution
% of wair to the Doppler velocity.
%
% Simon de Szoeke October 2009

visc_dyn=air_dyn_visc(T);
g=9.8;
rhow=1e3;
c=2*g*rhow/(9*visc_dyn);

SigX2=(2*log(DopplerV)-log(DopplerVWidth))/24;
SigX=sqrt(SigX2);
if nargout>1
    ModalRadius=sqrt(DopplerVfall.*exp(-14*SigX2)/c);
    if nargin>2 && nargout>2
        N=Z.*(2*ModalRadius).^(-6).*exp(-18*SigX2);
    end
end