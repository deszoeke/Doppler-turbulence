% examples of phase functions for different size parameters x
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

x=100;
m=1.4;
u=cosd([0:0.001:10 170.001:.001:180])';
