function y=ave_large2small(timesmall,timelarge,x)
% y=ave_large2small(timesmall,timelarge,x)
% Resamples long-time-interval data to short time intervals, copying the long-interval data
% to the short intervals as many times as sampled within the long interval.
% Timelarge must be a member of timesmall.
%
% Simon de Szoeke (c) 2009

% descandant of ave10to5.m

y=NaN+zeros(length(timesmall),size(x,2));
for i=1:length(timelarge)-1
    y(timesmall>=timelarge(i) & timesmall<timelarge(i+1))=x(i);
end
y(timesmall>=timelarge(end))=x(end); % last semiinfinite condition
return