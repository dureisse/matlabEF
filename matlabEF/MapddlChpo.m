function [mapddl1] = MapddlChpo(chpo1,nmail1,numer1,ListComp1);
% DUREISSEIX David    L.M.T. STRUCTURES et SYSTEMES  le 29 / 07 / 2002

% A partir d'un champ par point (chpo1,nmail1)
% d'une numerotation de nbno1 noeuds numer1 et d'une liste
% de nbco1 composantes ListComp1,
% construit une matrice de mapping mapddl1(nbno1,nbco1)
% mapddl1(no1,co1) = numero d'apparition dans la liste de ddl
% si 0, ddl inexistant

nbno1 = length(numer1);
nbco1 = length(ListComp1);
mapddl1 = zeros(nbno1,nbco1);

% On repere d'abord ou sont les ddl presents
% """""""""""""""""
% Boucle sur les sous-zones
nbzone1 = length(chpo1);
for zo1 = 1:nbzone1
  chpoe1 = chpo1{zo1};
  nmaile1 = nmail1{zo1};
  topo1 = nmaile1.MAIL;
  [junk,ia,ib] = intersect(nmaile1.MAIL,numer1);
%  junk=nmaile1.MAIL(ia)=numer1(ib)
  list_node1 = sort(ib);
% Boucle sur les composantes
  nbcomp1 = length(chpoe1);
  for comp1 = 1:nbcomp1
    [junk,ia,ib] = intersect(chpoe1{comp1}.COMP,ListComp1);
%   junkjunk=chpoe1{comp1}.COMP(ia)=ListComp1(ib)
    mapddl1(list_node1,ib) = 1;
  end
end

% On numerote ensuite
% """""""""""""""""""
ni = 0;
for no1 = 1:nbno1
  ddl1 = find(mapddl1(no1,:));
  nddl1 = length(ddl1);
  mapddl1(no1,ddl1) = mapddl1(no1,ddl1) + [ni:ni+nddl1-1];
  ni = ni + nddl1;
%  disp([int2str(no1) ...
%	' noeud ' int2str(numer1(no1)) ' nb ddl ' int2str(nddl1)])
end
