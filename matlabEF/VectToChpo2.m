function [chpo1,nmail1] = VectToChpo2(U1,numer1,mapddl1,ListDdl1, ...
                                      varargin)
% Disassemble a vector into a nodal field
%
% DUREISSEIX David    L.M.T. STRUCTURES et SYSTEMES  le 29 / 07 / 2002
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONATCTS  le 10 / 04 / 2004
%   Optional list of units
%
% Inputs
%  U1(nbddl1,1)			Vector
%  numer1(nbnot1)		List of nodes
%  mapddl1(nbnot1,nbcomp1)	Assembling mapping matrix
%  ListDdl1{nbcomp1}		Array of dof names
% Optional inputs
%  ListUnit1{nbcomp1}		Array of unit names
%                               (if not provided, '' is assumed)
% Outputs
%  chpo1,nmail1			Nodal field structures
%
% Comments
%   A partir d'un vecteur U1,
%   construit le champ par point (chpo1,nmail1)
%   de maniere "full" c'est a dire que tous les noeuds sont dans le
%   nuage de point sous-tendant nmail1, qu'il n'y a qu'une seule
%   sous-zone, et qu'il y a toutes les composantes.
%   Les valeurs inexistantes sont extrapolees a 0.
%   Ce format est celui attendu par AVS.

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
[nbnot1,nbcomp1] = size(mapddl1);

if ((nbcomp1 ~= length(ListDdl1)) || (nbcomp1 ~= length(ListUnit1)))
  nbcomp1
  length(ListDdl1)
  length(ListUnit1)
  error('unconsistent data')
end
if nbnot1 ~= length(numer1)
  nbnot1
  length(numer1)
  error('unconsistent data (bis)')
end

% Boucle sur les composantes
clear chpoel1;
for comp1 = 1:nbcomp1
  list_node1 = find(mapddl1(:,comp1));
  list_ddl1 = mapddl1(list_node1,comp1);
  xval1 = zeros(nbnot1,1);
  xval1(list_node1,1) = U1(list_ddl1,1);
  chpoel1{comp1} = struct('COMP',ListDdl1{comp1}, ...
                          'UNIT',ListUnit1{comp1}, ...
                          'XVAL',xval1);
end

% Une seule sous-zone
zo1 = 1;
clear chpo1 nmail1;

chpo1{zo1} = chpoel1;

% list_type_C3M1{4}='POI1'
nmail1{zo1} = struct('TYPE',list_type_C3M1{4},'MAIL',numer1');
