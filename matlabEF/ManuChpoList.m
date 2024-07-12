function [chpo1] = ManuChpoList(nmail1,ListComp1,ListUnit1,ListVal1)
% DUREISSEIX David  L.M.G.C. SYSTEMES MULTICONTACTS  le 13 / 04 / 2003          
%
% Construction d'un champ par point constant manuellement
% Idem ManuChpo mais syntaxe differente : ici on donne des listes

nvalue1 = length(ListComp1);
if (nvalue1 ~= length(ListUnit1) | nvalue1 ~= length(ListVal1))
  ListComp1
  ListUnit1
  ListVal1
  error('bad number of informations')
end

% Number of elements
nbzone0 = TestMeshType(nmail1,'POI1');
if (nbzone0 == 0)
  error('Mesh should be POI1')
end

% Loop on zones
clear chpo1;
for zo1 = 1:nbzone0
  nbpt1 = size(nmail1{zo1}.MAIL,1);
% Fill in
  clear chpoe1;
  for c1=1:nvalue1
    nomi   = ListComp1{c1};
    uniti  = ListUnit1{c1};
    valuei = ListVal1(c1);
    chpoe1{c1} = struct('COMP',nomi,'UNIT',uniti, ...
                        'XVAL',valuei*ones(nbpt1,1));
  end
  chpo1{zo1} = chpoe1;
end
