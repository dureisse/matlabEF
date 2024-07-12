function [mapddlDual1,mapddlPrimal1] = MapddlRigi2(modl1,mail1, ...
                                              numerd1,listCompDual1, ...
                                              numerp1,listCompPrimal1)
% DUREISSEIX David    L.M.T. STRUCTURES et SYSTEMES  le 02 / 08 / 2002
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 16 / 11 / 2002
%   une petite simplification : avec modl1, plus besoin de la rigidite
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 08 / 10 / 2007
%   Amelioration des performances

% A partir d'un modele (modl1,mail1)
% d'une numerotation duale de nbnod1 noeuds numerd1,
% d'une numerotation primale de nbnop1 noeuds numerp1,
% d'une liste duale et primale de nbcod1 et nbcop1 composantes
% listCompDual1 et listCompPrimal1,
% construit deux matrices de mapping mapddlDual1(nbnod1,nbcod1) et
% mapddlPrimal1(nbnop1,nbcop1).
% mapddlxxx1(no1,co1) = numero d'apparition dans la liste de ddl
% si 0, ddl inexistant
nbnod1 = length(numerd1);
nbnop1 = length(numerp1);
nbcod1 = length(listCompDual1);
nbcop1 = length(listCompPrimal1);
mapddlDual1 = zeros(nbnod1,nbcod1);
mapddlPrimal1 = zeros(nbnop1,nbcop1);


% Pour le dual (les lignes de la future matrice)
% ============

% On repere d'abord ou sont les ddl presents
% """""""""""""""""
% Boucle sur les sous-zones
nbzone1 = length(modl1);
for zo1 = 1:nbzone1
  modle1 = modl1{zo1};
  maile1 = mail1{zo1};

% On ne prend que les ddl de NDDD qui existent dans listCompDual1
  masq_pos2d = findoccur(modle1.DDLD(modle1.NDDD),listCompDual1);
% On ne prend que les ddl de NDDP qui existent dans listCompPrimal1
  masq_pos2p = findoccur(modle1.DDLP(modle1.NDDP),listCompPrimal1);
% On inverse une fois pour toute les listes
  numerd1_inv = InverseList(numerd1,max(numerd1));
  numerp1_inv = InverseList(numerp1,max(numerp1));

% Boucle sur les elements
  nbel1 = size(maile1.MAIL,1);
  for el1 = 1:nbel1
    topo1 = maile1.MAIL(el1,:);

%   On ne prend que les noeuds de NNOD qui existent dans numerd1
%%DD 07/10/08     masq_pos1d = findoccur(topo1(modle1.NNOD),numerd1);
    ii = max(topo1(modle1.NNOD)); if ii > length(numerd1_inv); numerd1_inv(ii) = 0; end
    masq_pos1d = numerd1_inv(topo1(modle1.NNOD));

%%DD 07/10/08 %   On ne prend que les ddl de NDDD qui existent dans listCompDual1
%%DD 07/10/08     masq_pos2d = findoccur(modle1.DDLD(modle1.NDDD),listCompDual1);

%   On ne prend que l'intersection des deux pour le dual
    masq3 = masq_pos1d & masq_pos2d;
    ld1 = find(masq3);
    pos_in_numerd1 = masq_pos1d(ld1);
    pos_in_compd1  = masq_pos2d(ld1);
    for i = 1:length(ld1)
      mapddlDual1(pos_in_numerd1(i),pos_in_compd1(i)) = 1;
    end

%   On ne prend que les noeuds de NNOP qui existent dans numerp1
%%DD 07/10/08     masq_pos1p = findoccur(topo1(modle1.NNOP),numerp1);
    ii = max(topo1(modle1.NNOP)); if ii > length(numerp1_inv); numerp1_inv(ii) = 0; end
    masq_pos1p = numerp1_inv(topo1(modle1.NNOP));

%%DD 07/10/08 %   On ne prend que les ddl de NDDP qui existent dans listCompPrimal1
%%DD 07/10/08     masq_pos2p = findoccur(modle1.DDLP(modle1.NDDP),listCompPrimal1);

%   On ne prend que l'intersection des deux pour le primal
    masq3 = masq_pos1p & masq_pos2p;
    lp1 = find(masq3);
    pos_in_numerp1 = masq_pos1p(lp1);
    pos_in_compp1  = masq_pos2p(lp1);
    for i = 1:length(lp1)
      mapddlPrimal1(pos_in_numerp1(i),pos_in_compp1(i)) = 1;
    end

  end

end

% On numerote ensuite
% """""""""""""""""""
ni = 0;
for no1 = 1:nbnod1
  ddl1 = find(mapddlDual1(no1,:));
  nddl1 = length(ddl1);
  mapddlDual1(no1,ddl1) = mapddlDual1(no1,ddl1) + [ni:ni+nddl1-1];
  ni = ni + nddl1;
%  disp([int2str(no1) ...
%	' noeud ' int2str(numerd1(no1)) ' nb ddl ' int2str(nddl1)])
end

ni = 0;
for no1 = 1:nbnop1
  ddl1 = find(mapddlPrimal1(no1,:));
  nddl1 = length(ddl1);
  mapddlPrimal1(no1,ddl1) = mapddlPrimal1(no1,ddl1) + [ni:ni+nddl1-1];
  ni = ni + nddl1;
%  disp([int2str(no1) ...
%	' noeud ' int2str(numerp1(no1)) ' nb ddl ' int2str(nddl1)])
end
