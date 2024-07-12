function [nbelt1] = Nelements2(mail1,varargin);
% Number of elements of a mesh (total or per zone)
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 30 / 12 / 2002
%
% Inputs
%   mail1		Mesh
% Optional inputs
%   listzone1(nzo1)	List of zones where to find the number of elt.
% Outputs
%   nbelt1		Number of elements in Mesh
% Comments
%   if listzone1 is not provided, nbelt1 is the total number of elts,
%   else, nbelt1(nzo1) is the list of element numbers in each zone.
%

nbzone1 = length(mail1);

narg = nargin-1;
if (narg == 0)
  listzone1 = [1:nbzone1];
elseif (narg == 1)
  listzone1 = varargin{1};
else
  narg
  error('Bar number of arguments')
end

clear nbelt1;
for izo1 = 1:length(listzone1)
  zo1 = listzone1(izo1);
  nbelt1(izo1) = size(mail1{zo1}.MAIL,1);
end

if (narg == 0)
  nbelt1 = sum(nbelt1);
end
