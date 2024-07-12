function [T1] = ChpoToChamnoOperator(numer2,mail2,numerptg);
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 22 / 03 / 2003
%
% Construit la matrice de l'operateur qui
% transforme un champ par point scalaire (1 composante) en un champ par
% elements defini au noeuds.
%
% Entrees
%   numer2(nbno2)	Numerotation locale des noeuds pour un vecteur
%                       representant le champ par point
%   mail2		Maillage cible
%   numerptg(Nptg1)     Numerotation locale des points support (i.e. les
%                       noeuds des elements)
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

  topo2 = reshape(topo2',1,nbel*nbno);
  topo2 = numer_inv(topo2);

  ptg2 = [ind1+1:ind1+nbno*nbel];
  ptg2 = numerptg_inv(ptg2);

  for i = 1:length(topo2)
    T1(ptg2(i),topo2(i)) = 1.;
  end
end

clear numer_inv numerptg_inv;
