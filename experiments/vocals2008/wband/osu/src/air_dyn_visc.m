function mu=air_dyn_visc(T)
% mu=air_dyn_visc(T)
% Compute the dynamic viscosity of air [Pa s] from the temperature [K]
% using Sutherland's formula.
%
% Simon de Szoeke 2009

T0=291.15; % K reference temperature
C=120;     % K Sutherland's constant for air
mu0=18.27e-6; % Pa s reference dynamic viscosity Pa s = kg m/s2/m2 s = kg/m/s

mu=mu0*(T0+C)./(T+C).*sqrt((T/T0).^3);

% multiply by density to get kinematic viscosity