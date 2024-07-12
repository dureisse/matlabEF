function [T1] = ChpoToChamOperator(numer2,mail2,numerptg,intg2);
% Operator that transform a nodal vector to an element vector
%
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 22 / 03 / 2003
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 04 / 07 / 2005
%
% Construit la matrice de l'operateur qui
% transforme un champ par point scalaire (1 composante) en un champ par
% elements.
%
% Entrees
%   numer2(nbno2)	Numerotation locale des noeuds pour un vecteur
%                       representant le champ par point
%   mail2		Maillage cible
%   numerptg(Nptg1)     Numerotation locale des points integration
%   intg2		Segment d'integration
%
% Sorties
%   T1(Nptg1,nbno2)	Matrice de l'operateur
%
% (voir ChpoToChamno)

% Number of target points
Nptg1 = length(numerptg); 
% Number of source points
nbno2 = length(numer2);
% Inverse renumbering
numer_inv = InverseList(numer2,max(numer2));
numerptg_inv = InverseList(numerptg,max(numerptg));

T1 = sparse(Nptg1,nbno2);

% Boucle sur les zones de mail2
ind1 = 0;
nbzone2 = length(mail2);
for zo2 = 1:nbzone2

  maile2 = mail2{zo2};
  topo2 = maile2.MAIL;
  [nbel,nbno] = size(topo2);

  intge2 = intg2{zo2};
  phi2   = intge2.PHI; % phi2(nbno,nbpgt)

  for el2 = 1:nbel
    topoe = topo2(el2,:);
    topoe = numer_inv(topoe);

    ptge = [ind1+(nbno*(el2-1))+1:ind1+(nbno*(el2-1))+nbno];
    ptge = numerptg_inv(ptge);

    T1(ptge,topoe) = phi2';
  end

  clear topo2 maile2;
end

clear numer_inv numerptg_inv;
