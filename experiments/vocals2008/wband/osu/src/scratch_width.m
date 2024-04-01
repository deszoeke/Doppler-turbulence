% Assuming Stokes v~r^2 fall velocity, find lognormal dropsize distribution
% width sigx from Doppler width/velocity with a lookup table.
sigx_stokes=0:.05:2;
ratio_stokes=2*exp(sigx.^2).*sqrt(cosh(sigx.^2).*sinh(sigx.^2));

% solution to find sigx
% interpolate sigx from Doppler width/velocity ratio using exact formula
sigxi=interp1(ratio_stokes,sigx_stokes,width./vel); % Stokes
sigxi=sqrt(log(1+width.^2/(vel+b/a).^2)); % linear velocity regime Fairall and 
% might need to compute for Gossard (1990) regime...hope not bc it's a small reduction of fall speed compared to linear

% approximations to the Stokes vel. width
apprx=2*exp(sigx.^2).*sigx;
apprx2=2*sigx.*(1+sigx.^2);
apprx3=2*sigx;

plot(sigx_stokes,ratio_stokes,'b.-')
hold on
plot(sigx, apprx,'r')
plot(sigx,apprx2,'m')
plot(sigx,apprx3,'k')
set(gca,'ylim',[0 5])

% Assume parameters of the lognormal distribution change only very slowly.
% Turbulent velocities influence the measured Doppler velocity fast,
% and should be ~0 in the long term (~10 min) average.