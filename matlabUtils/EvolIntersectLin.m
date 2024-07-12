function [x,y] = EvolIntersectLin(x1,y1,x2,y2);
% Evolutions linear intersection
%
% Inputs
%   x1(l1)	values of abscissa of first evolution
%   y1(l1)	values of first evolution
%   x2(l2)	values of abscissa of second evolution
%   y2(l2)	values of second evolution
%  Outputs
%   x(l)	values of abscissa of intersections
%   y(l)	values of intersections
%
% Linear interpolation is used
% Evolutions should truly intersect (i.e. change of sign on the
% difference of values)
% x1 and x2 should be increasing

if ~all(sort(x1)==x1) | ~all(sort(x2)==x2)
  error('not increasing')
end

lfin = unique([x1,x2]);
ii = find(lfin>=max(min(x1),min(x2)));
lfin = lfin(ii);
ii = find(lfin<=min(max(x1),max(x2)));
lfin = lfin(ii);

yfin1 = interp1(x1,y1,lfin);
yfin2 = interp1(x2,y2,lfin);

dy = yfin2 - yfin1;
ds = sign(dy) + 1;
ds(find(ds>0.5)) = 1;
dds = circshift(ds,[0,-1]);
dds = ds(1:end-1) - dds(1:end-1);

li2 = find(dds~=0);

x = ((lfin(li2) .* dy(li2+1)) - (lfin(li2+1) .* dy(li2))) ./ ...
    (dy(li2+1) - dy(li2));
y = interp1(x1,y1,x);
