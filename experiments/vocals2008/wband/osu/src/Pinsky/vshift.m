function y = vshift(x, kshift)
% shift vertically (dim 2) to align cloud tops

isf = isfinite(x);
sz = size(x);
nv = sz(2);

y = NaN(sz);
for i = 1:sz(1) % each time
    % shift dim 2 so top is at top index, preserving order
    %y(i,:) = circshift(x(i,:), [0, nv-ktop(i)]);
    y(i,:) = circshift(x(i,:), [0, kshift(i)]);
end
end
