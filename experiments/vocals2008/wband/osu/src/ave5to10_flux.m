function y=ave5to10(jdy5,jdy10,x)
% y=ave5to10(jdy5,jdy10,x)
% Average 5-min data x beginning at jdy5 to 10-minute intervals
% beginning at jdy10. Time is the leading dimension of x. y(i,:) is the
% average of x( jdy10(i)<=jdy5<jdy10(i+1) ,:).
% Warning: Differences between jdy5 and jdy10 will
% confuse the result of this function! jdy10 should be a member of jdy5.

% SPdeS
if length(jdy5)~=length(x)
    warning('AVE5TO10: length(jdy5)~=length(x)')
end

y=NaN+zeros(length(jdy10),size(x,2));
%isfin=isfinite(x);
for i=1:length(jdy10)-1
    ii=(jdy5>=jdy10(i) & jdy5<jdy10(i+1)); % [...)
    y(i,:)=nanmean(x(ii,:));
end
y(end,:)=nanmean(x(jdy5>=jdy10(end),:)); % semiinfinite condition [...+inf)
return
