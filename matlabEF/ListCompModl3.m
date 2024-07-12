function [listCompDual1,listCompPrimal1] = ListCompModl3 ...
	  (modl1,list_zones);
% DUREISSEIX David    L.M.T. STRUCTURES et SYSTEMES  le 13 / 11 / 2002
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 15 / 11 / 2002
%  modification des modeles version 1.3 passe a version 1.4

% Renvoie les noms de composantes primales et duales d'un modele
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
clear listCompPrimal1 listCompDual1;
  listCompPrimal1 = []; 
  listCompDual1   = [];

narg = nargin; % Nombre d'arguments
if (narg == 1)
  list_zones = [1:nbzone1];
elseif (narg == 2)
else
  error('Bar number of arguments')
end

ni = 0;
% Boucle sur les sous-zones du modele
for ii = 1:length(list_zones)
  zo1 = list_zones(ii);
  modle1 = modl1{zo1};

% Boucle sur les noms de composantes primales presentes dans le modele
  list_compe1 = modle1.COMP(modle1.NCOP);
  nbcomp1 = length(list_compe1);
  for i = 1:nbcomp1
    name1 = list_compe1{i};
	if ~isempty(listCompPrimal1)
    if ~ismember(name1,listCompPrimal1)
      ni = ni + 1;
      listCompPrimal1{ni} = name1;
    end
	else
      ni = ni + 1;
      listCompPrimal1{ni} = name1;
	end
  end
end

ni = 0;
% Boucle sur les sous-zones du modele
for ii = 1:length(list_zones)
  zo1 = list_zones(ii);
  modle1 = modl1{zo1};
% Boucle sur les noms de composantes duales presentes dans le modele
  list_compe1 = modle1.COMD(modle1.NCOD);
  nbcomp1 = length(list_compe1);
  for i = 1:nbcomp1
    name1 = list_compe1{i};
	if ~isempty(listCompDual1)
    if ~ismember(name1,listCompDual1)
      ni = ni + 1;
      listCompDual1{ni} = name1;
    end
	else
      ni = ni + 1;
      listCompDual1{ni} = name1;
	end
  end
end
