function P = Mie_phasefcn(m, x, u)
% Computation of Mie phase function for given 
% complex refractive-index ratio m=m'+im" 
% and size parameter x=k0*a, where k0= wave number in sphere 
% medium, a=sphere radius, using  Mie Scattering functions
% S1 and S2  and u=cos(scattering angle),
% where k0=vacuum wave number, a=sphere radius;
% s. p. 111-114, Bohren and Huffman (1983) BEWI:TDD122
% Simon de Szoeke adapted from C. Mätzler, May 2002

nmax=round(2+x+4*x^(1/3));
nang=length(u);
ab=Mie_ab(m,x);
an=ab(1,:);
bn=ab(2,:);

[p,t]=Mie_pt_vect(u,nmax); %p,t (length(u),nmax)

n=(1:nmax);
n2=(2*n+1)./(n.*(n+1));
pin=ones(nang,1)*n2.*p;
tin=ones(nang,1)*n2.*t;

S1=(an*pin'+bn*tin');
S2=(an*tin'+bn*pin');

% scattering
q=2*(n+1)*(abs(an).^2+abs(bn).^2)';
%qsca=2*q/x^2;

% phase function
P=(abs(S1).^2+abs(S2).^2)'/q;
return

% Plot examples:
m=1.4;
u=cosd([0:.5:15 16:1:44 45:5:180])'
clf
for ix=-2:2;
    x=10^ix;
    P=Mie_phasefcn(m,x,u);
    h(ix+3)=polar(acos(u),P/max(P));
    hold on
end
set(h(2),'color','g')
set(h(3),'color','k')
set(h(4),'color','r')
set(h(5),'color','m')
legend('x=0.01','x=0.1','x=1','x=10','x=100')
title 'Mie phase functions'
