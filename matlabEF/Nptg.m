function [nbptg1] = Nptg(mail1,intg1,varargin)
% Number of integration points (total or per zone)
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 02 / 08 / 2002
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 12 / 04 / 2004
%  Ajout argument optionnel
%
% Inputs
%   mail1               Mesh
%   intg1		Integration information
% Optional inputs
%   listzone1(nzo1)     List of zones where to find the number of ptg
% Outputs
%   nbptg1              Number of integration points in Mesh
% Comments
%   if listzone1 is not provided, nbptg1 is the total number of ptg,
%   else, nbptg1(nzo1) is the list of ptg numbers in each zone.

nbzone1 = length(mail1);
if (nbzone1 ~= length(intg1))
  nbzone1
  length(intg1)
  error('Unconsistent data')
end

narg = nargin-2;
if (narg == 0)
  listzone1 = [1:nbzone1];
elseif (narg == 1)
  listzone1 = varargin{1};
else
  narg
  error('Bar number of arguments')
end

clear nbptg1;
for izo1 = 1:length(listzone1)
  zo1 = listzone1(izo1);
  nbptg1(zo1) = size(mail1{zo1}.MAIL,1) * size(intg1{zo1}.WEIGHT,2);
end

if (narg == 0)
  nbptg1 = sum(nbptg1);
end
