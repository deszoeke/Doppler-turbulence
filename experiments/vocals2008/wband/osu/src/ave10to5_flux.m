function y=ave10to5(jdy5,jdy10,x)
% y=ave10to5(jdy5,jdy10,x)
% Resamples 10-min data to 5-min intervals, copying the low rate data as
% many times as sampled within the low rate interval.
% Warning: jdy10 must be a member of jdy5.

% SPdeS
y=NaN+zeros(length(jdy5),size(x,2));
for i=1:length(jdy10)-1
    y(jdy5>=jdy10(i) & jdy5<jdy10(i+1))=x(i);
end
y(jdy5>=jdy10(end))=x(end); % last semiinfinite condition
return