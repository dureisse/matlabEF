function [Lchpo1,nmail1] = LvectToLchpo2(LU1, ...
                                         numer1,mapddl1,ListDdl1, ...
                                         varargin)
% Disassemble a set of vectors into a list of nodal fields
%
% DUREISSEIX David    L.M.T. STRUCTURES et SYSTEMES  le 27 / 10 / 2003
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONATCTS  le 10 / 04 / 2004
%   Optional list of units
%
% Entrees
%   LU1(nddl,n)			Liste de n vecteurs
%   numer1(nbnot1)		Numerotation des noeuds
%   mapddl1(nbnot1,nbcomp1)	Matrice de mapping pour assemblage
%   ListDdl1{nbcomp1}		Liste des noms de ddl
% Entrees optionnelles
%   ListUnit1{nbcomp1}		Liste des unites
% Sorties
%   Lchpo1{n}			Array de la liste de champs par point
%   nmail1			Maillage POI1 sous tendant
%
% Comments
%   Meme routine que VectToChpo2, mais travaillant sur des listes :
%   a partir d'une liste de vecteurs LU1, construit la liste de
%   champs par point (Lchpo1,nmail1)
%   de maniere full c'est a dire que tous les noeuds sont dans le
%   nuage de point sous-tendant nmail1, qu'il n'y a qu'une seule
%   sous-zone, et qu'il y a toutes les composantes.
%   Les valeurs inexistantes sont extrapolees a 0.
%   Ce format est celui attendu par AVS.
%
%   Tous les champs par point d'une liste ont la meme structure.

nin1 = nargin-4;
if nin1 == 0
  ListUnit1 = repmat([{''}],length(ListDdl1));
elseif nin1 == 1
  ListUnit1 = varargin{1};
else
  nin1
  error('Wrong number of optional arguments')
end

GlobalVar;
[nddl,n] = size(LU1);
[nbnot1,nbcomp1] = size(mapddl1);

if ((nbcomp1 ~= length(ListDdl1)) || (nbcomp1 ~= length(ListUnit1)))
  nbcomp1
  length(ListDdl1)
  length(ListUnit1)
  error('no consistent lengths')
end
if nbnot1 ~= length(numer1)
  nbnot1
  length(numer1)
  error('no consistent lengths bis')
end

% Une seule sous-zone
zo1 = 1;
% list_type_C3M1{4}='POI1'
clear nmail1;
nmail1{zo1} = struct('TYPE',list_type_C3M1{4},'MAIL',numer1');


% Boucle sur les champs
clear Lchpo1;
for i = 1:n

% Boucle sur les composantes
  clear chpoel1;
  for comp1 = 1:nbcomp1
    list_node1 = find(mapddl1(:,comp1));
    list_ddl1 = mapddl1(list_node1,comp1);
    xval1 = zeros(nbnot1,1);
    xval1(list_node1,1) = LU1(list_ddl1,i);
    chpoel1{comp1} = struct('COMP',ListDdl1{comp1}, ...
			    'UNIT',ListUnit1{comp1}, ...
			    'XVAL',xval1);
    clear list_node1 list_ddl1 xval1;
  end

  clear chpo1;
  chpo1{zo1} = chpoel1;
  Lchpo1{i} = chpo1;
  clear chpoel1 chpo1;
end
