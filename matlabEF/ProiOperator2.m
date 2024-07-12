function [T1,notfound] = ProiOperator2(mail1,numerptg,xcoor2,xcoor1,xcrit1);
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 17 / 03 / 2003
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACT  le 16 / 03 / 2004
%   Ajout de la liste notfound
%
% Construit un operateur de projection discretise entre
% vecteurs de champs scalaires definis sur des
% maillages ou des points d'integration differents.
%
% Il permet de passer d'un vecteur representant un champ par
% elements defini aux noeuds a un vecteur representant
% un champ par point en des noeuds particuliers,
% en utilisant les fonctions de forme de l'element geometrique
% pour interpoler.
% (voir Proi)
%
% Entrees
%   mail1		Maillage support du champ par element aux noeuds
%   numerptg(Nptg1,1)	Numerotation locale des points support (i.e. les
%                       noeuds dans les elements) 
%   xcoor2(nbno2,idim)	Coordonnees des points cibles
%   xcoor1(:,idim)	Pile des noeuds
%   xcrit1  : critere relatif de proximite (sera multiplie par la taille
%             de l'element courant pour avoir le critere absolu)
%
% Sorties
%   T1(nbno2,Nptg1)	Operateur de projection
%   notfound		Liste des points cibles non trouves
%
% Remarque : les noeuds du vecteur representant un champ par point
% sont ranges dans l'ordre de xcoor2.


% Number of integration points
Nptg1 = length(numerptg);
% Number of target points and dimension
[nbno2,idim] = size(xcoor2);

notfound = [];

% Transformation matrix between vectors  U1 = T1 * S1
% """""""""""""""""""""""""""""""""""""""""""""""""""

% reperage des noeuds deja trouves
dejatrouv2 = zeros(1,nbno2);
toto = zeros(1,nbno2);

T1 = sparse(nbno2,Nptg1);

% Loop on elements of mail1
iptg = 0;
nbzone1 = length(mail1);
for zo1 = 1:nbzone1
  maile1 = mail1{zo1};
  type1  = maile1.TYPE;
  Xcorel1 = EF_CoorRefNod(type1);
  topo1 = maile1.MAIL;
  [nbel,nbno] = size(topo1);
  for el1 = 1:nbel
    if nbel > 10
      disp(['  element ' int2str(el1) ' / ' int2str(nbel) ...
            ' in zone ' int2str(zo1)])
    end

%   Find not-yet-found candidates (inside the box)
    xcorel1 = xcoor1(topo1(el1,:),:);
    max1 = max(xcorel1);
    min1 = min(xcorel1);
    eps1 = xcrit1 * max(max1 - min1);
    test1 = (xcoor2 > repmat(min1-eps1,nbno2,1)) & ...
            (xcoor2 < repmat(max1+eps1,nbno2,1));
    lpot1 = find(all(test1,2)' & ~dejatrouv2);

%   Find candidates that are really inside the element
    for ipot1 = 1:length(lpot1)
      locno2 = lpot1(ipot1);
      ino2 = locno2;
      xcor2 = xcoor2(ino2,:);

      [idans,xcor1,shp1,xcor2a] = EF_CoorRef3( ...
                         type1,Xcorel1,xcorel1,xcor2,eps1,1.2*xcrit1);
      if idans
        ListPtg = [iptg+1:iptg+nbno];
        ListLocPtg = findoccur(ListPtg,numerptg);
        ind = find(ListLocPtg);
        ListLocPtg = ListLocPtg(ind);
        T1(locno2,ListLocPtg) = shp1(:,ind);
        dejatrouv2(locno2) = 1;
        toto(locno2) = norm(xcor2a-xcor2);
      end
      clear xcor1 shp1 xcor2a xcor2 ino2 locno2;
    end

    iptg = iptg + nbno;
    clear lpot1 test1 xcorel1;
  end
  clear maile1 topo1;
end
n_reste = length(find(dejatrouv2 == 0));
if n_reste
  disp(['ProiOperator: Warning ' int2str(n_reste) ' nodes not found'])
  toto(n_reste)
  notfound = [notfound find(dejatrouv2 == 0)];
end
