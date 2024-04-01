function d = dmn(x, indx)
[a, n] = accu(x, indx);
d = a./n;
d(n<5) = NaN; % QC the means a bit
end

function [a, n] = accu(x, indx)
ii = isfinite(x);
xi = x;
xi(~ii) = 0;
a = zeros(length(unique(indx))+1, size(x,2));
n = zeros(size(a));
for i = 1:length(indx)
    a(indx(i),:) = a(indx(i),:) + xi(i,:);
    n(indx(i),:) = n(indx(i),:) + ii(i,:);
end
end
