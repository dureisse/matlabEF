function [chml2] = ExtrChml2(chml1,ListComp1,varargin)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 28 / 12 / 2002
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 23 / 09 / 2002
%   Possibilite de changer le nom de la composante en arg. optionnel

% Extrait le champ par element constant chml2 a partir du
% champ par element constant chml1 constitue des seules
% composantes presentes dans ListComp1(ncomp1)
% Si ListComp2(ncomp1) est fourni, les noms des ncomp1 composantes
% sont changes en consequence.

narg = nargin; % Nombre d'arguments
if (narg == 3)
  ListComp2 = varargin{1};
else
  ListComp2 = ListComp1;
end

nbcomp2 = length(ListComp1);
if (nbcomp2 ~= length(ListComp2))
  ListComp1
  ListComp2
  error('bad number of new component names')
end

List1 = ListCompChml(chml1);
ListZone1 = findoccur(ListComp1,List1);
if length(find(ListZone1)) ~= nbcomp2
  ListComp1
  List1
  error('Not all the components where found')
end

clear chml2;
for zo2 = 1:length(ListZone1)
  zo1 = ListZone1(zo2);
  chml2{zo2} = chml1{zo1};
  chml2{zo2}.COMP = ListComp2{zo2};
end
