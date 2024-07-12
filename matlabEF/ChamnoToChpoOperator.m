function [S,D] = ChamnoToChpoOperator(modl2,mail2,intg2, ...
                                      numer2,numer1,xcoor1);
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 13 / 04 / 2003
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 07 / 04 / 2005
%  Correction de quelques bugs
%
% Construit l'operateur discretise LI = S * D qui transforme un
% champ par elements defini au noeuds en un champ par point
% Champs a 1 seule composante
% Operation de lissage - ponderation par la surface des elements
% (voir ChamnoToChpo)
%
% Entrees
%   modl2		Modele (sert a trouver les noeuds de la
%                               transformation)
%   mail2		Maillage du champ par elements defini aux noeuds
%   intg2		Segment d'integration aux noeuds
%   numer2(Nptg2)	Numerotation des points d'integration (aux noeuds)
%   numer1(Nno1)	Numerotation des noeuds du champ par point
%   xcoor1(nbnot,idim)	Pile des noeuds
%
% Sorties
%  S(Nno1,Nptg2)	Operateur de sommation aux noeuds
%  D(Nptg2,Nptg2)	Matrice de ponderation
%


if (1 == 0)
error('Cette routine est erronee... voir ChamnoToChpo')
% On construit l'operateur discret de sommation aux noeuds S
% """"""""""""""""""""""""""""""""""""""""""""""""""""""""""
Nno1 = length(numer1);
Ncomp2 = 1;
Nptg2 = length(numer2);
S = sparse(Nno1,Nptg2);

numer_inv2 = InverseList(numer2,max(numer2));          
numer_inv1 = InverseList(numer1,max(numer1));          

ind2 = 0;
nbzone2 = length(mail2);
for zo2 = 1:nbzone2
  maile2 = mail2{zo2};
  topo2 = maile2.MAIL;
  [nbel2,nbno2] = size(topo2);
  list_ptg2 = [ind2+1:ind2+nbel2*nbno2];
  list_ptg2 = numer_inv2(list_ptg2);
  list_no1 = reshape(topo2',1,nbel2*nbno2);
  list_no1 = numer_inv1(list_no1);
  for i = 1:nbel2*nbno2
    S(list_no1(i),list_ptg2(i)) = 1.;
  end

  ind2 = ind2 + nbel2*nbno2;
  clear list_no1 list_ptg2 topo2 maile2;
end

end


% On construit l'operateur de duplication aux ptg
% """""""""""""""""""""""""""""""""""""""""""""""
T1 = ChpoToChamnoOperator(numer1,mail2,numer2);

% Remarque : l'operateur discret de sommation aux noeuds est T1'
% """"""""


% On construit le chamno discret de valeur aire de l'element A2
% """""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Nptg2 = Nptg(mail2,intg2);
A2 = zeros(Nptg2,1);

mode1 = 'COPL';
idim = size(xcoor1,2);
[modl2,intg2a] = ModlIntg13(mail2,'ELASTIQUE','RIGIDITE',mode1,idim,intg2);
mode1 = '';
[r1,modlr1] = IntegrOperator2(modl2,mail2,intg2a,xcoor1,mode1);

numer_inv2 = InverseList(numer2,max(numer2));
ind2 = 0;
nbzone2 = length(mail2);
for zo2 = 1:nbzone2
  maile2 = mail2{zo2};
  topo2 = maile2.MAIL;
  [nbel2,nbno2] = size(topo2);

  re1 = r1{zo2};
  xval1 = re1.XVAL;

  for el2 = 1:nbel2
    aire2 = sum(sum(xval1(:,:,el2)));
    list_ptg2 = [ind2+1:ind2+nbno2];
    list_ptg2 = numer_inv2(list_ptg2);
    A2(list_ptg2,1) = aire2;

    ind2 = ind2 + nbno2;
    clear list_ptg2;
  end

  clear xval1 re1 topo2 maile2;
end
clear r1 modlr1;

% On construit la matrice de ponderation D
% """"""""""""""""""""""""""""""""""""""""
%D = diag(((T1*S*A2) .^ (-1.)) .* A2);
D = diag(((T1*(T1'*A2)) .^ (-1)) .* A2);

% On a l'operateur de lissage LI
% """"""""""""""""""""""""""""""
%LI = S * D;
%LI = T1' * D;

S = T1';
D = sparse(D);
