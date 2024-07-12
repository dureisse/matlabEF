function [max1] = MaxiChml(chml1,varargin)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 29 / 12 / 2002
%
% max1 = Maxichml(chml1)
% max1 = Maxichml(chml1,ListComp1)
% Trouve la valeur maximale d'un champ par element constant chml1
% pour les composantes dont les noms sont dans le deuxieme
% argument s'il existe (sinon, partout).

narg = nargin; % Nombre d'arguments
if (narg == 1)
  ListComp1 = ListCompChml(chml1);
elseif (narg == 2)
  ListComp1 = varargin{1};
else
  narg
  error('Bad number of arguments')
end

max1 = [];
nbcomp1 = length(chml1);
for i1 = 1:nbcomp1
  chmle1 = chml1{i1};
  l1 = findoccur([{chmle1.COMP}],ListComp1);
  if l1
    max1 = max([max1 max(max(chmle1.XVAL))]);
  end
  clear chmle1;
end
