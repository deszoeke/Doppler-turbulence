function y=ave_small2large(timesmall,timelarge,x)
% y=ave_small2large(timesmall,timelarge,x)
%
% Average data x on small time intervals beginning at timesmall to longer intervals
% beginning at timelarge. Time is the leading dimension of x. y(i,:) is the
% average of x( timelarge(i)<=timesmall<timelarge(i+1) ,:).
%
% Timelarge should be a member of timesmall. Small differences
% between timesmall and timelarge cause uneven results from the
% inequality deciding the averaging interval.
%
% Simon de Szoeke (c) 2009

% Descendant of ave5to10.m

if length(timesmall)~=length(x)
    warning('AVE_SMaLL2LARGE: length(timesmall)~=length(x)')
end

y=NaN+zeros(length(timelarge),size(x,2));
%isfin=isfinite(x);
for i=1:length(timelarge)-1
    ii=(timesmall>=timelarge(i) & timesmall<timelarge(i+1)); % [...)
    y(i,:)=nanmean(x(ii,:));
end
y(end,:)=nanmean(x(timesmall>=timelarge(end),:)); % semiinfinite condition [...+inf)
return
