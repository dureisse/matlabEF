function [vec3] = PscalVect(vec1,numer1,map1,list1, ...
                            vec2,numer2,map2,list2, ...
                            listps1,listps2,numer3)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 03 / 10 / 2003
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 04 / 10 / 2003
%   Cas particulier d'un champ scalaire
%
% Effectue le produit scalaire de deux champs representes sous forme
% de vecteurs. Les composantes multipliees entre elles sont celles
% donnees dans les listes listps1 et listps2 (nom a nom).
% listps1 et listps2 doivent donc etre de meme longueur.
% Le cas particulier ou listps1 (ou listps2) n'a qu'une composante est
% celui ou le champ correspondant est scalaire ; on multiplie alors
% le deuxieme champ par la valeur scalaire du premier.
%
% Inputs
%   vec1(n1,1)		Premier vecteur
%   numer1(nbno1)	Numerotation des noeuds pour le mapping
%   map1(nbno1,nbco1)	Matrice de mapping
%   list1(nbco1)	Liste des noms de composantes ou ddl
%   vec2(n2,1)		Deuxieme vecteur
%   numer2(nbno2)	Numerotation des noeuds pour le mapping
%   map2(nbno2,nbco2)	Matrice de mapping
%   list2(nbco2)	Liste des noms de composantes ou ddl
%   listps1(nbco)	Liste des composantes "primales"
%   listps2(nbco)	Liste des composantes "duales"
%   numer3(nbno3)	Numerotation des noeuds pour construire le
%                       vecteur resultat
% Outputs
%  vec3(nbno3,1)	Vecteur resultat
%
% Pour le mapping, la valeur correspondant a la composante de nom
% list1(co1) pour le noeud numer1(no1) est vec1(map1(no1,co1)).
% Comme le vecteur resultat vec3 est scalaire, il n'a qu'une
% composante qui n'a pas a etre precisee. Son assemblage se fait
% donc uniquement dans l'ordre de numer3. Son mapping serait
% map3 = [1:nbno3]'
% Dans le cas d'un champ scalaire, le mapping est le meme que celui
% du champ non scalaire.
% Si des noeuds ne sont pas dans numer3, il ne sont pas assembles;
% si des noeuds de numer3 n'existent pas dans vec1 ou vec2, on met
% la valeur nulle.


nbno3 = length(numer3);

if length(listps1) == length(listps2)

% Produit scalaire
  vec3 = zeros(nbno3,1);

% Composantes existantes
  lco1 = findoccur(listps1,list1); ico1 = (lco1 ~= 0);
  lco2 = findoccur(listps2,list2); ico2 = (lco2 ~= 0);
  ico = find(ico1 & ico2);
  lco1 = lco1(ico);
  lco2 = lco2(ico);

  for no3 = 1:nbno3
    no1 = findoccur(numer3(no3),numer1);
    no2 = findoccur(numer3(no3),numer2);
    if ((no1 ~= 0) && (no2 ~= 0))
%     Noeud existant
      vec3(no3,1) = vec1(map1(no1,lco1),1)' * vec2(map2(no2,lco2),1);
    end
  end

elseif length(listps1) == 1

% Produit exterieur par un scalaire
  nbco = length(listps2);
  vec3 = zeros(nbno3*nbco,1);

% Composantes existantes
  lco1 = findoccur(listps1,list1); ico1 = (lco1 ~= 0);
  if isempty(ico1)
    return;
  end
  lco2 = findoccur(listps2,list2); ico2 = (lco2 ~= 0);
  ico = find(ico2);
  lco2 = lco2(ico);

  for no3 = 1:nbno3
    no1 = findoccur(numer3(no3),numer1);
    no2 = findoccur(numer3(no3),numer2);
    if (no1 ~= 0) & (no2 ~= 0)
%     Noeud existant
      vec3(map2(no2,lco2),1) = vec1(map1(no1,lco1),1)' * ...
           vec2(map2(no2,lco2),1);
    end
  end

elseif length(listps2) == 1
  error('Pas encore prevu, desole... retournez les arguments')

else
  listps1
  listps2
  error('Bad number of components')
end
