function [mapcom1,mapddl1] = MapcomB2(modl1,mail1,intg1, ...
                                      numerptg1,listCom1, ...
                                      numer1,listDdl1, ...
                                      opti1)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 19 / 11 / 2002
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 08 / 10 / 2007
%  Amelioration des performances

% A partir d'un modele (modl1,mail1,intg1)
% d'une numerotation points d'integration numerptg1 (numerotation globale),
% d'une liste de noms de composantes listCom1,
% d'une numerotation de noeuds numer1, et
% d'une liste de noms de ddl listDdl1,
% construit deux matrices de mapping mapcom1 et mapddl1.
% (voir MapddlRigi2)
% opti1 est une option ('DUAL' ou 'PRIMAL') qui indique 
% si on prend les composantes et ddl primaux ou duaux.

nbno1 = length(numer1);
nbptg1 = length(numerptg1);
nbcom1 = length(listCom1);
nbddl1 = length(listDdl1);
mapcom1 = zeros(nbptg1,nbcom1);
mapddl1 = zeros(nbno1,nbddl1);

% Pour les ddl
% ============
%
% On repere d'abord ou sont les ddl et les composantes presents
% """""""""""""""""
ind = 0; % for point numbering

% Boucle sur les sous-zones
nbzone1 = length(modl1);
for zo1 = 1:nbzone1
  modle1 = modl1{zo1};
  maile1 = mail1{zo1};
  intge1 = intg1{zo1};

  nbptg = length(intge1.WEIGHT);                                                

  if strcmp(opti1,'PRIMAL')
    namec1 = modle1.COMP(modle1.NCOP);
    named1 = modle1.DDLP(modle1.NDDP);
    nnox = modle1.NNOP;
  elseif strcmp(opti1,'DUAL')
    namec1 = modle1.COMD(modle1.NCOD);
    named1 = modle1.DDLD(modle1.NDDD);
    nnox = modle1.NNOD;
  else
    opti1
    erreur('option unknown')
  end
% On ne prend que les composantes de NCO? qui existent dans listCom1
  name1f    = repmat(namec1,1,nbptg);
  masq_pos4 = findoccur(name1f,listCom1);
% On ne prend que les ddl de NDD? qui existent dans listDdl1
  masq_pos2 = findoccur(named1,listDdl1);
% On inverse une fois pour toute les listes
  numer1_inv    = InverseList(numer1,max(numer1));
  numerptg1_inv = InverseList(numerptg1,max(numerptg1));

% Boucle sur les elements
  nbel1 = size(maile1.MAIL,1);
  for el1 = 1:nbel1
    topo1 = maile1.MAIL(el1,:);
    topo2 = [ind+1:ind+nbptg]; 
    ind = ind + nbptg;

%   On ne prend que les noeuds de NNO? qui existent dans numer1
%%DD 07/10/08    masq_pos1 = findoccur(topo1(nnox),numer1);
    ii = max(topo1(modle1.NNOD)); if ii > length(numer1_inv); numer1_inv(ii) = 0; end
    masq_pos1 = numer1_inv(topo1(modle1.NNOD));

%   On ne prend que l'intersection avec les composantes existantes
    masq3 = masq_pos1 & masq_pos2;
    ld1 = find(masq3);
    pos_in_numer1 = masq_pos1(ld1);
    pos_in_comp1  = masq_pos2(ld1);
    for i = 1:length(ld1)
      mapddl1(pos_in_numer1(i),pos_in_comp1(i)) = 1;
    end

%   On ne prend que les points qui existent dans numerptg1
    topo2f = reshape(repmat(topo2,nbcom1,1),1,nbcom1*nbptg);
%%DD 07/10/08    masq_pos1 = findoccur(topo2f,numerptg1);
    ii = max(topo2f); if ii > length(numerptg1_inv); numerptg1_inv(ii) = 0; end
    masq_pos1 = numerptg1_inv(topo2f);

%   On ne prend que l'intersection avec les composantes existantes              
    masq3 = masq_pos1 & masq_pos4;
    lp1 = find(masq3);
    pos_in_numerptg1 = masq_pos1(lp1);
    pos_in_com1  = masq_pos4(lp1);
    for i = 1:length(lp1)
      mapcom1(pos_in_numerptg1(i),pos_in_com1(i)) = 1;
    end

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
ni = 0;
for no1 = 1:nbptg1
  ddl1 = find(mapcom1(no1,:));
  nddl1 = length(ddl1);
  mapcom1(no1,ddl1) = mapcom1(no1,ddl1) + [ni:ni+nddl1-1];
  ni = ni + nddl1;
%  disp([int2str(no1) ...
%	' noeud ' int2str(numerptg1(no1)) ' nb ddl ' int2str(nddl1)])
end
