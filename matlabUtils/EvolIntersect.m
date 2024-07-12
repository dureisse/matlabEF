function [x,y] = EvolIntersect(x1,y1,x2,y2,varargin);
% Evolutions intersection
%
% Inputs
%   x1(l1)	values of abscissa of first evolution
%   y1(l1)	values of first evolution
%   x2(l2)	values of abscissa of second evolution
%   y2(l2)	values of second evolution
%   option      optional word, i.e. 'lin' or 'semilogy' or 'loglog'
%               for interpolation (default value is 'lin')
%  Outputs
%   x(l)	values of abscissa of intersections
%   y(l)	values of intersections
%
% Linear interpolation is used
% Evolutions should truly intersect (i.e. change of sign on the
% difference of values)
% x1 and x2 should be increasing

nin1 = nargin-4;
if (nin1 == 0)
  option = 'lin';
else
  option = varargin{1};
end

if strcmp(option,'lin')
 [x,y] = EvolIntersectLin(x1,y1,x2,y2);
elseif strcmp(option,'semilogy')
 [x,y] = EvolIntersectLin(x1,log(y1),x2,log(y2));
 y = exp(y);
elseif strcmp(option,'loglog')
 [x,y] = EvolIntersectLin(log(x1),log(y1),log(x2),log(y2));
 y = exp(y);
 x = exp(x);
else
   option
  error('bad option')
end
