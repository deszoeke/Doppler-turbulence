% SPdeS
% Demonstrate that a bimodal distribution can have different widths depending on
% how the peaks are chosen.

v=-2:0.1:8;
sig1=0.5;
v1=1;
sig2=1;
v2=3;
g1=2*exp(-(v-v1).^2/(2*sig1.^2));
g2=1*exp(-(v-v2).^2/(2*sig2.^2));
plot(v,g1,v,g2,v,g1+g2)

[gmin,timin]=min(g1(v>v1 & v<v2)+g2(v>v1 & v<v2)); % local min.
imin=timin+find(v==v1);

% first peak, left of local min
ii=v<v(imin);
[g1max,i1max]=max(ii.*g1+g2);
v1peak=v(i1max);
width1 =sqrt(sum(((v(ii)-v1peak).^2).*(g1(ii)+g2(ii)))/sum(g1(ii)+g2(ii)));
% compare to sig1

% second peak, right of local min
ii=v>v(imin);
[g2max,i2max]=max(ii.*g1+g2);
v2peak=v(i2max);
width2 =sqrt(sum(((v(ii)-v2peak).^2).*(g1(ii)+g2(ii)))/sum(g1(ii)+g2(ii)));
% compare to sig2

v12=sum(v.*(g1+g2))/sum(g1+g2);
width12=sqrt(sum((v.^2).*(g1+g2))/sum(g1+g2)-v12^2);