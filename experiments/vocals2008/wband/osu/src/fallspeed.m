function [output, a, b]=fallspeed(input,sw,alt)
% speed=fallspeed(radius)   or   radius=fallspeed(speed,-1);
% 
% Compute the droplet fall speed (m/s) for a drop radius (m)
% or radius given the fall speed, if a negative 2nd argument is given.
% using the function in Gossard et al. 1990. Updated for 1.5 km and 3 K, e.g.
% stratus clouds (VOCALS 2008). Uses air_dyn_visc(T) to compute the viscosity
% for the Stokes drop fall speed. Nominal altitude [km] may be specified
% after the inverse switch, e.g. speed=fallspeed(radius,1,1.5); default
% altitude is 1.5 km. Optional arguments 2 and 3 are the linear fit parameters
% a and b for r=aV+b between the Stokes and the Mason regime (V: [0.4 2.75] m/s)
% 
% VOCALS 2008 :: 2009-09-25 :: Simon de Szoeke

if ~exist('alt');
    alt=1.5; % km
elseif alt<0
    disp('Warning: altitude <0 km')
end
% parameters for 1.5 km, 276.45 K
g=9.8; % m/s2
rho0=1.19;
gamma=9.8; % K/km
T0=291.15; % K (18 C)
T=T0-gamma*alt;

rhow=1e3;
rho=rho0*exp(-alt/8);
rhofac=exp(-0.4*alt/8); % (rho/rho0)^.4
visc_dyn=air_dyn_visc(T); % 0.0175

% a=1.1e-4;  % s  linear best fit to other curves
% b=1.25e-5; % m
% Solve for linear relation connecting Vf=[0.4 2.75] m/s
sp=[0.4 2.75];
ra(1)=sqrt(9*visc_dyn.*sp(1)/(2*g*rhow)); % m
ra(2)=0.5e-3*(-1.667)*log((9.65-sp(2).*rhofac)/10.3); %m
a=(ra(2)-ra(1))/(sp(2)-sp(1));
b=ra(1)-a*sp(1);
% optionally return a and b

if ~exist('sw','var')
    sw=1;
end

if sw<0 % invert: solve for radius from speed
    speed=abs(input);
    i1=speed<sp(1);
    i2=speed>=sp(1) & speed<sp(2);
    i3=speed>=sp(2);
    if ~isempty(i1)    % Stokes
        radius(i1)=sqrt(9*visc_dyn.*speed(i1)/(2*g*rhow)); % m
    end
    if ~isempty(i2)    % linear relation
        radius(i2)=a*speed(i2)+b; % m
    end
    if ~isempty(i3)
        radius(i3)=0.5e-3*(-1.667)*log((9.65-speed(i3).*rhofac)/10.3); %m
    end
    output=radius;
else
    radius=abs(input);
    i1=radius<sqrt(9*visc_dyn.*splim(1)/(2*g*rhow));% speed<0.4
    i2=radius>=sqrt(9*visc_dyn.*splim(1)/(2*g*rhow)) & ... % speed>=0.4 & speed<2.75
        radius<-1.667*log((9.65-splim(2).*rhofac)/10.3);
    i3=radius>=-1.667*log((9.65-splim(2).*rhofac)/10.3);
    i3=speed>=2.75;
    if ~isempty(i1)    % Stokes
        speed(i1)=2*r.^2*g*rhow/(9*visc_dyn); % m/s
    end
    if ~isempty(i2)    % linear blend
        speed(i2)=(radius-b)/a; % m/s
    end
    if ~isempty(i3)
        speed(i3)=(9.65-10.3*exp(2*radius/-1.667e-3))/rhofac; % m/s
    end
    output=speed;
end
return

if false
% scratch:
V1=[[0 1 3 6]*1e-3 0.01:0.01:0.4]'; % m/s, Stokes regime
V2=(0.4:0.05:2.8)'; % linear
V3=[1:0.2:10.4]'; % Mason fit

D1=sqrt(18*visc_dyn*1e3.*V1/g);                    % Stokes
D2=0.22*V2 + 0.022;                          % linear - visual match
DGossard=0.2*V2;
DFrisch=(1.2e-4*V2+1e-5)*2e3;
D3=-1.667*log((9.65-V3.*rhofac)/10.3);     % rain
% new way, radius in meters...
R1=sqrt(9*visc_dyn.*V1/(2*g*rhow)); % m
R2=a*V2+b;
R3=0.5e-3*(-1.667)*log((9.65-V3.*rhofac)/10.3); % unchanged

plot(R1,V1,'-','linewidth',2)
hold on
%plot(DGossard,V2,'.-','color',[0 0.5 0])
%plot(DFrisch,V2,'.-','color',[0.5 0 0.5])
plot(R2,V2,'k-','linewidth',2)
plot(R3,V3,'-','color',[1 0 0],'linewidth',2)
set(gca,'yscale','log','xscale','log','xlim',[1e-6 2e-3],'ylim',[1e-3 10],'fontsize',16)
xlabel 'radius (m)'
ylabel 'fall speed (m/s)'

end