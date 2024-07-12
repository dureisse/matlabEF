function [S1] = ChamToVect3(cham1,mail1,intg1,numer1,mapComp1,listComp1)
% DUREISSEIX David    L.M.T. STRUCTURES et SYSTEMES  le 03 / 08 / 2002
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 29 / 08 / 2002

% Place un champ (non constant) par element (cham1,mail1,intg1)
% dans un vecteur S1
% selon le mapping mapComp1 (voir MapddlCham)
% associe a numer1 et listComp1.

% L'ordre des composantes est celui de listComp1, 
% La numerotation locale des points d'integration
% est celle du maillage mail1, on prend uniquement les
% points de numer1.

% DE GROSSES SIMILARITES AVEC ChpoToVect3...
% LES REGROUPER ??
% mail1 NE SERT PAS

% point = noeud pour Chpo
%       = point d'integration pour Cham

nbddl1 = max(max(mapComp1));
S1 = zeros(nbddl1,1);

ind = 0;

% Boucle sur les sous-zones du champ par element
nbzone1 = length(cham1);
for zo1 = 1:nbzone1
  chamel1 = cham1{zo1};
%%  maile1  = mail1{zo1};
  intge1  = intg1{zo1};
  nbpt1   = length(intge1.WEIGHT);
% Boucle sur les composantes du champ par element
%%  nbcomp1 = size(chamel1,2);
  nbcomp1 = length(chamel1);
  for i = 1:nbcomp1
    xval1 = chamel1{i}.XVAL;
    [nbelt,nbvalt] = size(xval1);
    xval1t = reshape(xval1',nbelt*nbvalt,1);
    nbval  = nbvalt / nbpt1;
    if (nbval ~= 1)
      error('multiple value component not implemented')
    end

%   Component in listComp1
    j = findoccur({chamel1{i}.COMP},listComp1);

%   Points (d'integration) in cham1, in numer1
%%DD 27/01/03    list_pt1 = find(mapComp1(:,i));
%%DD 25/08/03    list_pt1 = find(mapComp1(:,j));
    j = j(find(j)); list_pt1 = find(mapComp1(:,j));
    list_pt1 = find(mapComp1(:,j));
    list_pt3 = numer1(list_pt1);  % points in numer1
%%DD 27/01/03    list_pt2 = [ind+1:ind+nbpt1]; % points in subzone of mail1
%%DD 27/01/03    ind = ind + nbpt1;
    list_pt2 = [ind+1:ind+nbpt1*nbelt]; % points in subzone of mail1
    [junk,ia,ib] = intersect(list_pt3,list_pt2);
%   junk=list_pt3(ia)=list_pt2(ib)
%   list_pt3(ia)=numer1(list_pt1(ia))

    list_ddl1  = mapComp1(list_pt1(ia),j);
%%DD 27/01/03    S1(list_ddl1,1) = S1(list_ddl1,1) + xval1(ib,1);
    S1(list_ddl1,1) = S1(list_ddl1,1) + xval1t(ib,1);
  end
  ind = ind + nbpt1*nbelt;
end
