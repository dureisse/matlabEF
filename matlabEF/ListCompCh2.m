function [listComp1,listUnit1] = ListCompCh2(ch1,list_zones);
% DUREISSEIX David    L.M.T. STRUCTURES et SYSTEMES  le 31 / 07 / 2002

% Renvoie les composantes d'un champ (par element ou par point)
% et leur unite (voir AVS)
% Si list_zones n'existe pas, sur toutes les zones
% sinon dans celles de la liste list_zones

nbzone1 = size(ch1,2);
ni = 0;
clear listComp1 listUnit1;
  listComp1 = []; listUnit1 = [];

narg = nargin; % Nombre d'arguments
if (narg == 1)
  list_zones = [1:nbzone1];
elseif (narg == 2)
else
  error('Bar number of arguments')
end

% Boucle sur les sous-zones du champ
for ii = 1:length(list_zones)
  zo1 = list_zones(ii);
  chel1 = ch1{zo1};
% Boucle sur les composantes du champ
  nbcomp1 = size(chel1,2);
  for i = 1:nbcomp1
    name1 = chel1{i}.COMP;
    if isempty(listComp1)
      ni = ni + 1;
      listComp1{ni} = chel1{i}.COMP;
      listUnit1{ni} = chel1{i}.UNIT;
    elseif ~ismember(name1,listComp1)
      ni = ni + 1;
      listComp1{ni} = chel1{i}.COMP;
      listUnit1{ni} = chel1{i}.UNIT;
    end
  end
end
