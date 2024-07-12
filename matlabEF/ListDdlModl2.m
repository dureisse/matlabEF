function [listDdlDual1,listDdlPrimal1] = ListDdlModl2 ...
	  (modl1,list_zones);
% DUREISSEIX David    L.M.T. STRUCTURES et SYSTEMES  le 13 / 11 / 2002

% Renvoie les noms de ddl primaux et duaux d'un modele
% Si list_zones n'existe pas, sur toutes les zones
% sinon dans celles de la liste list_zones

narg = nargin; % Nombre d'arguments
if (narg == 1)
  list_zones = [1:length(modl1)];
elseif (narg == 2)
else
  error('Bar number of arguments')
end



nbzone1 = length(modl1);
clear listDdlPrimal1 listDdlDual1;
  listDdlPrimal1 = []; 
  listDdlDual1 = [];

narg = nargin; % Nombre d'arguments
if (narg == 1)
  list_zones = [1:nbzone1];
elseif (narg == 2)
else
  error('Bar number of arguments')
end

ni = 0;
% Boucle sur les sous-zones de la rigidite
for ii = 1:length(list_zones)
  zo1 = list_zones(ii);
  modle1 = modl1{zo1};
% Boucle sur les noms de ddl primaux presents dans le modele
  list_ddle1 = modle1.DDLP(unique(modle1.NDDP));
  nbddl1 = length(list_ddle1);
  for i = 1:nbddl1
    name1 = list_ddle1{i};
    if ~isempty(listDdlPrimal1)
    if ~ismember(name1,listDdlPrimal1)
      ni = ni + 1;
      listDdlPrimal1{ni} = name1;
    end
    else
      ni = ni + 1;
      listDdlPrimal1{ni} = name1;
    end
  end
end

ni = 0;
% Boucle sur les sous-zones de la rigidite
for ii = 1:length(list_zones)
  zo1 = list_zones(ii);
  modle1 = modl1{zo1};
% Boucle sur les noms de ddl duaux presents dans le modele
  list_ddle1 = modle1.DDLD(unique(modle1.NDDD));
  nbddl1 = length(list_ddle1);
  for i = 1:nbddl1
    name1 = list_ddle1{i};
    if ~isempty(listDdlDual1)
    if ~ismember(name1,listDdlDual1)
      ni = ni + 1;
      listDdlDual1{ni} = name1;
    end
    else
      ni = ni + 1;
      listDdlDual1{ni} = name1;
    end
  end
end
