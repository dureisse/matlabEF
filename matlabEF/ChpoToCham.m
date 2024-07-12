function [cham2] = ChpoToCham(chpo1,nmail1,mail2,intg2);
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 02 / 04 / 2003
%
% Transforme un champ par point (chpo1,nmail1) en un champ par
% elements (cham2,mail2,intg2)
%
% Remarque : les composante non definies partout sur mail2 sont
% extrapolees par 0.

% On fait un pseudo assemblage sur les noeuds de mail2
nmail2 = ChangeMesh2(mail2,'POI1'); % Une seule zone
numer2 = nmail2{1}.MAIL';
numer_inv2 =  InverseList(numer2,max(numer2));
clear nmail2;
[listComp1,listUnit1] = ListCompChpo2(chpo1);
mapddl1 = MapddlChpo(chpo1,nmail1,numer2,listComp1);
C1 = ChpoToVect3(chpo1,nmail1,numer2,mapddl1,listComp1);


nbcomp1 = length(listComp1);
nbzone2 = length(mail2);
clear cham2;

% Boucle sur les zones de mail2
for zo2 = 1:nbzone2
  clear chame2;
  maile2 = mail2{zo2};
  intge2 = intg2{zo2};
  topo2 = maile2.MAIL;
  nod2 = numer_inv2(topo2);
  phi2 = intge2.PHI;

% Boucle sur les composantes de chpo1
  for i = 1:nbcomp1
    lmap1 = mapddl1(:,i)';

%   Valeurs aux noeuds
%subtil
    xval = zeros(size(nod2));
    xval(:) = C1(lmap1(nod2));

%   Valeurs aux ptg
%   On suppose autant de fct de formes que de noeuds et dans l'ordre
    xval = xval * phi2;

    chame2{i} = struct('COMP',listComp1{i},'UNIT',listUnit1{i}, ...
                       'XVAL',xval);
    clear lmap1 xval;
  end

  cham2{zo2} = chame2;
  clear nod2 topo2 maile2 chame2;
end

clear C1 mapddl1 listComp1 listUnit1 numer_inv2 numer2;
