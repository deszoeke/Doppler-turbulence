function [SigX, ModalRadius, N]=drizzle_lognorm_param(DopplerVfall, DopplerVfallWidth, Z)
% [SigX, ModalRadius, N]=drizzle_lognorm_param(DopplerVfall, DopplerVfallWidth, Z);
%  ln(m)   m                                        m/s            m/s          m^3
%
% Find the parameters SigX. ModalRadius, and N for lognormally distributed
% drizzle drops,
%
% n_drizzle(x) = N/(SigX*sqrt(2*pi))*exp[-(x-x0)^2/2*SigX^2]
%
% with x=ln(r), x0=ln(ModalRadius), from a single realization Doppler drop
% fall speed and Doppler fall speed width. A linear fall-speed radius relation for
% drizzle is adopted r = a*DopplerVfall + b; a=1.2e-4 s, b=1.0e-5 m (Frisch et al.
% 1995). The Doppler radar velocity is actually the sum of the air w and the fall
% speed. The air velocity adjustment to the fall speed is not done by this function
% which expects FALL velocities and widths. Z and N are optional.
%
% Simon de Szoeke September 2009

a=1.2e-4; %s
b=1.0e-5; %m
SigX2=log( 1 + (DopplerVfallWidth./DopplerVfall + b/a).^2 ); % Eq. 14 Frisch 1995
SigX=sqrt(SigX2);
if nargout>1
    ModalRadius=(a*DopplerVfall+b).*exp(-13*SigX2/2);
    if nargin>2 && nargout>2
        N=Z.*(2*ModalRadius).^(-6).*exp(-18*SigX2);
    end
end