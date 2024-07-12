function [cham2] = ExtrCham2(cham1,ListComp1,varargin)
%   Extract a component of an element-based field
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 27 / 05 / 2004
%
% Extrait le champ par element cham2 a partir du
% champ par element cham1 constitue des seules
% composantes presentes dans ListComp1(ncomp1)
% Si ListComp2(ncomp1) est fourni, les noms des ncomp1 composantes
% sont changes en consequence.
%
% Voir ExtrChml2

narg = nargin; % Nombre d'arguments
if (narg == 3)
  ListComp2 = varargin{1};
else
  ListComp2 = ListComp1;
end

nbcomp2 = length(ListComp1);
if (nbcomp2 ~= length(ListComp2))
  ListComp1
  ListComp2
  error('bad number of new component names')
end

clear cham2;
nbzone1 = length(cham1);
for zo1 = 1:nbzone1
  chame1 = cham1{zo1};
  ncomp1 = length(chame1);
  chame2 = [];
  for i = 1:ncomp1
    j = findoccur({chame1{i}.COMP},ListComp1);
    if j
      ii = length(chame2);
      chame2{ii+1} = chame1{i};
      chame2{ii+1}.COMP = ListComp2{j};
    end
  end
  cham2{zo1} = chame2;
  clear chame1 chame2;
end
