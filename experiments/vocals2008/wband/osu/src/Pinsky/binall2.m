function indx = binall2( bin, x )
% indx = binall2(bin,x)
% BINALL2 : "Bin into bins"
% Returns index of bin into which the independent variable x
% falls. Bins are defined by semiopen intervals [bin(i) bin(i+1))
% Uses all indices of x so output indx has same size as x.
% NaNs and out-of-range values are mapped to indx = length(bin).
%
% Performs no operations, but can be used to bin-average data (see binavg).
%
% See also BIN2
% Simon de Szoeke, from binavg

% BIN2 uses only finite x
%isf=isfinite(x);
%xf=x(isf);
% But BINALL2 uses all x so size(indx) = size(x).

indx=floor(interp1(bin,1:length(bin), x, 'linear'));
%make a special bin at the end to catch data not in any bin
indx(indx<1 | ~isfinite(indx))=length(bin);
