function [mapcomp1] = MapCompCham(cham1,mail1,intg1,numer1,ListComp1);
% DUREISSEIX David    LMGC SYSTEMES MULTICONTACTS  le 29 / 08 / 2002

% A partir d'un champ par element (non constant) (cham1,mail1,intg1)
% d'une numerotation de nbpt1 points numer1 (points d'integration de
% la numerotation locale de mail1) et d'une liste
% de nbco1 composantes ListComp1,
% construit une matrice de mapping mapcomp1(nbpt1,nbco1)
% mapcomp(pt1,co1) = numero d'apparition dans la liste de composantes
% si 0, inexistant

nbpt1 = length(numer1);
nbco1 = length(ListComp1);
mapcomp1 = zeros(nbpt1,nbco1);


% On repere d'abord ou sont les composantes presentes
% """""""""""""""""
% Boucle sur les sous-zones
ind = 0;
nbzone1 = length(cham1);
for zo1 = 1:nbzone1
  chame1 = cham1{zo1};
  maile1 = mail1{zo1};
  intge1 = intg1{zo1};
%%%  topo1  = maile1.MAIL;
  nbptg1 = size(maile1.MAIL,1) * size(intge1.WEIGHT,2);
  topo1  = [ind+1:ind+nbptg1];
  ind = ind + nbptg1;
  [junk,ia,ib] = intersect(topo1,numer1);
%  junk=nmaile1.MAIL(ia)=numer1(ib)
  list_pt1 = sort(ib);
% Boucle sur les composantes
  nbcomp1 = length(chame1);
  for comp1 = 1:nbcomp1
    [junk,ia,ib] = intersect(chame1{comp1}.COMP,ListComp1);
%   junk=chame1{comp1}.COMP(ia)=ListComp1(ib)
    mapcomp1(list_pt1,ib) = 1;
  end
end

% On numerote ensuite
% """""""""""""""""""
ni = 0;
for pt1 = 1:nbpt1
  comp1  = find(mapcomp1(pt1,:));
  ncomp1 = length(comp1);
  mapcomp1(pt1,comp1) = mapcomp1(pt1,comp1) + [ni:ni+ncomp1-1];
  ni = ni + ncomp1;
%  disp([int2str(pt1) ...
%	' point ' int2str(numer1(pt1)) ' nb comp ' int2str(ncomp1)])
end
