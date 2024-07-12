function [chpo2,nmail2] = ExtrChpo(chpo1,nmail1,ListComp2,varargin)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 12 / 08 / 2003
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 24 / 12 / 2003
%  Changement de nom optionnel
%
% Extrait du champ par point chpo1, le champ par point
% chpo2 constitue des seules composantes presentes dans ListComp2(ncomp1) 
% En option, changement possible du nom de ces composantes
%
% Entrees
%   chpo1		Champ par point de depart
%   nmail1		Son maillage support
%   ListComp2(ncomp1)	Liste des composantes a extraire
% Entree optionnelle
%   ListComp3(ncomp1)	Liste des nouveaux noms
% Sorties
%   chpo2		Champ par point resultat
%   nmail2		Son maillage support
%
% Voir ExtrChml

if nargin == 3
  ListComp3 = ListComp2;
elseif nargin == 4
  ListComp3 = varargin{1};
  if length(ListComp3) ~= length(ListComp2)
    ListComp2
    ListComp3
    error('Bar list of names')
  end
else
  nargin
  error('Bar number of arguments')
end

clear chpo2 nmail2;
zo2 = 0;
% Loop on zones
nbzone1 = length(chpo1);
for zo1 = 1:nbzone1
  ListComp1 = ListCompChpo2(chpo1,zo1);
%%DD 03/12/25  list1 = findoccur(ListComp1,ListComp2);
  list1 = findoccur(ListComp2,ListComp1);
  list2 = find(list1);
  list1 = list1(list2);
  if ~isempty(list1)
%   If the component has been found, add a new zone and copy
    zo2 = zo2 + 1;
    chpoe1  = chpo1{zo1};
    nmaile1 = nmail1{zo1};
    clear chpoe2 nmaile2;
    for i = 1:length(list1)
      j1 = list1(i);
      j2 = list2(i);
      chpoe2{i} = chpoe1{j1};
      chpoe2{i}.COMP = ListComp3{j2}; 
    end
    nmaile2 = nmaile1;
    chpo2{zo2}  = chpoe2;
    nmail2{zo2} = nmaile2;
    clear chpoe1 chpoe2 nmaile1 nmaile2;
  end
  clear list1 ListComp1;
end
