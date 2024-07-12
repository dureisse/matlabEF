function [ListComp1] = ListCompChml(chml1)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 28 / 12 / 2002

% Extrait du champ par element constant chml1
% la liste des noms des composantes ListComp1(ncomp1)

nbcomp1 = length(chml1);
ListComp1 = [];

for zo1 = 1:nbcomp1
  ListComp1 = [ListComp1 {chml1{zo1}.COMP}];
end
