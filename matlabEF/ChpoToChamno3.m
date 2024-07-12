function [chamno2,intg1] = ChpoToChamno3(chpo1,nmail1,mail2);
% Transforms a nodal field to an element field defined at nodes
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 23 / 01 / 2003
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 22 / 03 / 2003
%   Ajout du segment d'integration intg1
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 20 / 01 / 2005
%   Traite les champs par point a plusieurs zones
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 26 / 02 / 2007
%   Correction d'un cas particulier pour le POI1
%
% Transforme un champ par point (chpo1,nmail1) en un champ par
% elements defini au noeuds (chamno2,mail2,intg2)
%
% Remarque : les composante non definies partout sur mail2 sont
% extrapolees par 0.
%
% TRES SIMILAIRE A ChpoToCham : PREVOIR DE LES UNIFIER PLUS QUE CA !
% PAR EXEMPLE DANS ChpoToCham : SI INTG2 NON FOURNI, ALORS
% ON FAIT COMME ChpoToChamno2 ET ON RETOURNE intg1
% PREVOIR AUSSI LE CAS ChpoToChml (il faudra peut etre un 
% SegmentIntgGravite)

% Pour permettre les champ par points discontinus d'une zone
% a l'autre, on boucle sur les zones de chpo1
nbzone1 = length(chpo1);
if nbzone1 > 1
  disp('  Treatment of possible discontinuous nodal fields')
end
[listComp1,listUnit1] = ListCompChpo2(chpo1);
nbcomp1 = length(listComp1);

% Pour faire un pseudo assemblage sur les noeuds de mail2
nmail2 = ChangeMesh2(mail2,'POI1'); % Une seule zone
numer2 = nmail2{1}.MAIL';
numer_inv2 =  InverseList(numer2,max(numer2));
clear nmail2;

% On cree le champ par elements vide
clear chamno2;
% Boucle sur les zones de mail2
nbzone2 = length(mail2);
for zo2 = 1:nbzone2
% Boucle sur les composantes de chpo1
  for i = 1:nbcomp1
    maile2 = mail2{zo2};
    topo2 = maile2.MAIL;
    nod2 = numer_inv2(topo2);
if (size(topo2,2)) == 1
% On corrige un effet parasite DD 26/02/2007
  nod2 = nod2';
end
    xval = zeros(size(nod2));
    chamnoe2{i} = struct('COMP',listComp1{i},'UNIT',listUnit1{i}, ...
                         'XVAL',xval);
    clear xval nod2 topo2 maile2;
  end
  chamno2{zo2} = chamnoe2;
  clear chamnoe2;
end


% Boucle sur les zones du champ par point source
% """"""""""""""""""""""""""""""""""""""""""""""
for zo1 = 1:nbzone1

% Assemblage de la zone uniquement sur tous les noeuds
% """"""""""""""""""""""""""""""""""""""""""""""""""""
  clear chpo0; chpo0{1} = chpo1{zo1};
  clear nmail0; nmail0{1} = nmail1{zo1};
  mapddl1 = MapddlChpo(chpo0,nmail0,numer2,listComp1);
  C1 = ChpoToVect3(chpo0,nmail0,numer2,mapddl1,listComp1);

% On somme le resultat dans le champ par element
% """"""""""""""""""""""""""""""""""""""""""""""
% Boucle sur les zones de mail2
  nbzone2 = length(mail2);
  for zo2 = 1:nbzone2
    maile2 = mail2{zo2};
    topo2 = maile2.MAIL;
    nod2 = numer_inv2(topo2);
if (size(topo2,2)) == 1
% On corrige un effet parasite DD 26/02/2007
  nod2 = nod2';
end
    chamnoe2 = chamno2{zo2};

%   Boucle sur les composantes de chpo1
    for i = 1:nbcomp1
      lmap1 = mapddl1(:,i)';
      xval = chamnoe2{i}.XVAL;
%      xval = xval + C1(lmap1(nod2));
toto = lmap1(nod2);
if (size(topo2,2)) == 1
% On corrige un effet parasite DD 26/02/2007
  toto = toto';
end

[ii,jj,vv] = find(toto);
for tt = 1:length(ii)
  xval(ii(tt),jj(tt)) = xval(ii(tt),jj(tt)) + C1(vv(tt));
end
clear ii jj tt toto;
      chamnoe2{i}.XVAL = xval;
      clear xval lmpa1;
    end
    chamno2{zo2} = chamnoe2;

    clear nod2 topo2 maile2;
  end

  clear C1 mapddl1 nmail0 chpo0;
end





% Pseudo integration segment associated to mail2
intg1 = SegmentIntgNo(mail2);

clear numer_inv2 numer2;
