% Use a Mie calculation to determine the Mie(more exact) to
% Rayleigh(more approx.) backscatter ratio at 94 GHz.

L=load('lhermitte_tableii.txt');
T=L(:,1); % temperature
M=L(:,2)+1i*L(:,3); % index of refraction for water at 94 GHz, BHMie expects m=n+ik n,k>0
K2=abs((M.^2-1)./(M.^2+2)).^2;
lambda=3.2e6;  % 3.2 mm in nm
k=2*pi/lambda;
r=10.^(-1:.0001:4)*0.5e3;            % exponentially distributed radius vector, nm
X=k*r;             % size parameter
nx=length(X);

nt=length(T);
    
if false
    qext=NaN+zeros(nx,nt);
    qsca=NaN+zeros(nx,nt);
    qabs=NaN+zeros(nx,nt);
    qb=NaN+zeros(nx,nt);
    asy=NaN+zeros(nx,nt);
    qratio=NaN+zeros(nx,nt);
    if exist('mie_efficiencies_94GHz.mat','file');
        load mie_efficiencies_94GHz
    else
        for it=1:4 % temperature
            for ix=1:nx
                % calculate with code from Maetzler's Mie.m:
                result=Mie(M(it),X(ix));
                qext(ix,it)=result(1);   % extinction efficiency
                qsca(ix,it)=result(2);   % scattering
                qabs(ix,it)=result(3);   % absorption
                qb(ix,it)=result(4);     % backscatter <<< use this one for reflectivity
                asy(ix,it)=result(5);    % asymetry parameter
                qratio(ix,it)=result(6); %
            end
        end
        save mie_efficiencies_94GHz T M X qext qsca qabs qb asy qratio
    end
end

% interpolate dielectric constants from Lhermitte to 9 C
t=9;
m=interp1(T,M,t);
k2=abs((m.^2-1)./(m.^2+2)).^2;

qext=NaN+zeros(nx,1);
qsca=NaN+zeros(nx,1);
qabs=NaN+zeros(nx,1);
qb=NaN+zeros(nx,1);
asy=NaN+zeros(nx,1);
qratio=NaN+zeros(nx,1);
if exist('mie_efficiencies_94GHz9C.mat','file');
    load mie_efficiencies_94GHz9C
else
    for ix=1:nx
        % calculate with code from Maetzler's Mie.m:
        result=Mie(m,X(ix));
        qext(ix)=result(1);   % extinction efficiency
        qsca(ix)=result(2);   % scattering
        qabs(ix)=result(3);   % absorption
        qb(ix)=result(4);     % backscatter <<< use this one for reflectivity
        asy(ix)=result(5);    % asymetry parameter
        qratio(ix)=result(6); %
    end
    save mie_efficiencies_94GHz9C m X qext qsca qabs qb asy qratio
end

loglog(2*r/1e3,qb,'b')
hold on
loglog(2*r/1e3,qext,'r')
set(gca,'fontsize',18)

gammaf=lambda^4*qb./(4*pi^4*k2*((2*r).^4))';
loglog(2*r/1e3,gammaf,'k')
xlabel('diameter 10^{-6} m')
ylabel('Mie backscatter, extinction effic. and Mie/Rayleigh backscatter ratio')
title('3.2 mm W-band radar scattering off liquid water spheres')
print -dpng 
% distributions with range of median volume radii R0
N=100; % irrelevant efficiencies and lidar ratio
R0=([5:1:90 100:10:900 1000:100:2000])'*1e3;  %10-2000 microns in nm

% lognormal distributions have median volume radius r0=1.444*rcent, the
% central (mode) radius rcent.
x=log(r);
sigx=0.35;
gammaprime=NaN+zeros(size(R0));
% loop over central radius
dr=[diff(r) 0];
const=lambda^4/(4*pi^4*k2);
for i=1:length(R0)
    r0=R0(i);
    x0=log(r0/1.444);
    n=N/(sigx*2*pi)*exp(-(x-x0).^2/(2*sigx.^2)); % lognormal radius spectrum
    % Mie backscatter coefficient yields Mie to Rayleigh backscatter ratio
    top=n.*(2*r).^2.*dr*qb; % summation in matrix multiplication inner product
    bot=n.*(2*r).^6*dr';
    gammaprime(i)=const*top/bot;
end

% gamma distributions
% (already transformed to be around median volume radius)
mu=[0 2 5 10]; % shape parameter, 2-10 preferred by Miles et al. 2000
nmu=length(mu);
fmu=6/3.67^4*(3.67+mu).^(mu+4)./gamma(mu+4); % ^(mu+4) corrected eqn 20 O'Connor (2005)
lwcg=NaN+zeros(size(R0,1),nmu);
ggprime=NaN+zeros(size(R0,1),nmu);
% loop over shape parameter mu and median volume radius R0
for imu=1:nmu; % slow!
    for i=1:length(R0)
        r0=R0(i);
        g=N*fmu(imu)*(r/r0).^mu(imu).*exp(-(3.67+mu(imu))*r/r0);
        top=g.*(2*r).^2.*dr*qb; % summation in matrix multiplication inner product
        bot=g.*(2*r).^6*dr';
        lwcg(i,imu)=4/3*pi*sum(g.*r.^3.*dr); % volume LWC
        ggprime(i,imu)=const*top/bot;
    end
end
semilogx(2*R0/1e3,gammaprime,'k.-')
hold on
semilogy(2*R0/1e3,ggprime)
axis([10 3e3 0 1.2])
set(gca,'fontsize',18)
xlabel('median volume diameter D_0')
ylabel('Mie-to-Rayleigh backscatter ratio \gamma''')
legend('lognormal','\mu=0','\mu=2','\mu=5','\mu=10')
title('94 GHz')
print('-dpng','Mie2Rayleigh_ratio.png')

% save Mie to Rayleigh ratios in each distribution in lookup table.
save gammaprime R0 mu gammaprime ggprime
