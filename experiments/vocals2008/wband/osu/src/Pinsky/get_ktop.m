function ktop = get_ktop(x)
sz = size(x);
isf = isfinite(x);
ktop = zeros(sz(1),1);
fii = find( any(isf,2) );
for j = 1:length(fii) % each time
     i = fii(j);
     ktop(i) = find(isf(i,:), 1,'last');
end
end
