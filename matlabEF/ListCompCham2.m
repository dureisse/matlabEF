function [listComp1,listUnit1] = ListCompCham2(cham1,list_zones);
% DUREISSEIX David    L.M.T. STRUCTURES et SYSTEMES  le 31 / 07 / 2002

% Renvoie les composantes d'un champ (par element ou par point)
% et leur unite (voir AVS)
% Si list_zones n'existe pas, sur toutes les zones
% sinon dans celles de la liste list_zones

narg = nargin; % Nombre d'arguments
if (narg == 1)
  [listComp1,listUnit1] = ListCompCh2(cham1);
elseif (narg == 2)
  [listComp1,listUnit1] = ListCompCh2(cham1,list_zones);
else
  error('Bar number of arguments')
end

